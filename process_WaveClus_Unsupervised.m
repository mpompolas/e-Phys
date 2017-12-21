function varargout = process_separate_RAW_per_electrode( varargin )

% @=============================================================================

% This process separates the initial raw signal to nChannels binary signals
% and performs spike sorting individually on each channel with the WaveClus
% spike-sorter. The spikes are clustered and assigned to individual
% neurons. The code ultimately produces a raw_elec(i)_spikes.mat
% for each electrode that can be used later for supervised spike-sorting.
% When all spikes on all electrodes have been clustered, all the spikes for
% each neuron is assigned to an events file in brainstorm format.

% The input is the number of samples that sequentially will be added to
% each file. The larger the number, the faster the procedure. Beware of
% memory issues for a large number.

% Konstantinos Nasiotis
% @=============================================================================




% Ways to improve this code:


% This was tested on .ns5 files (Blackrock - Ripple).
% These files are integers with int16 presision (uV) that are converted to
% (V). Other acquisition systems might not be in uV, or not be in int16
% precision.
% GENERALIZE



% Waveclus only reads txt or .mat files.
% I first create and append a binary file and later I convert it to a .mat.
% The mat files are inputs to the spike sorter which just detects where the
% spikes occured. These spikes need to be clustered after.


%     10 minutes recording at 30k Hz. 192 channels
% 1. Separating raw file to individual channels: 200 seconds
% 2. Converting binary files to .mat files: 230 seconds
% 3. Spike sorting: 520 seconds
% 4. Clustering: 430 seconds





% @=============================================================================
% This function is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2017 University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Francois Tadel, 2011-2014

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'WaveClus spike sorting';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Import recordings';
    sProcess.Index       = 1200;
    sProcess.Description = 'www.in.gr';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'raw'};
    sProcess.OutputTypes = {'raw'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.isSeparator = 1;
    sProcess.options.binsize.Comment = 'Samplesize for appending: ';
    sProcess.options.binsize.Type    = 'value';
    sProcess.options.binsize.Value   = {4, 'million samples', 1};
    sProcess.options.paral.Comment = 'Parallel processing';
    sProcess.options.paral.Type    = 'checkbox';
    sProcess.options.paral.Value   = 0;
    sProcess.options.make_plots.Comment = 'Create Images';
    sProcess.options.make_plots.Type    = 'checkbox';
    sProcess.options.make_plots.Value   = 0;
    % Channel name comment
    sProcess.options.make_plotshelp.Comment = '<I><FONT color="#777777">This saves images of the clustered spikes</FONT></I>';
    sProcess.options.make_plotshelp.Type    = 'label';
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    OutputFiles = {};
    [~, iStudy, ~] = bst_process('GetOutputStudy', sProcess, sInputs);
    sTargetStudy = bst_get('Study', iStudy);
    
    
    % Before doing anything, check if wave_clus exists is in the path
    if ~exist('wave_clus_font')
         bst_report('Error', sProcess, sInputs, 'This process requires the WaveClus spike-sorter.');
         return;
    end
    
    
    
    % Start computations   
    bst_progress('start', 'WaveClus spike-sorting', ['Demultiplexing raw file'], 0, length(sInputs));

    for i = 1:length(sInputs)
        bst_progress('start', 'WaveClus spike-sorting', 'Demultiplexing raw file...');
        tic
        DataFile = file_fullpath(sInputs(i).FileName); % E:\brainstorm_db\Playground\data\Monkey\@rawtest_LFP_EYE\data_0raw_test_LFP_EYE.mat
        ChannelFile = sInputs(i).ChannelFile;
        ChannelMat = in_bst_channel(ChannelFile); 
        
        % This is Francois's snippet for working on the temporary folder
        % Work in Brainstorm's temporary folder
        workDir = bst_fullfile(bst_get('BrainstormTmpDir'), ['Unsupervised_Spike_Sorting\' sInputs(i).Condition(5:end)]);
        % Make sure Matlab is not currently in the work directory
        curDir = pwd;
        if ~isempty(strfind(pwd, workDir))
            curDir = bst_fileparts(workDir);
            cd(curDir);
        end
        % Erase if it already exists
        if file_exist(workDir)            
            file_delete(workDir, 1, 3);
            pause(3)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % MARTIN
            
            if file_exist(workDir)
                STOP % THIS MEANS THE FOLDER WAS NOT DELETED BECAUSE A FILE IS STILL OPEN IN IT. THIS IS DANGEROUS CAUSE THAT FILE WILL BE APPENDED LATER IN THE CODE.
                     % GIVE A WARNING TO THE USER TO RESTART MATLAB or run
                     % the function again.
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end
        % Create empty work folder
        res = mkdir(workDir);
        if ~res
            bst_report('Error', sProcess, sInputs, ['Cannot create temporary directory: "' workDir '".']);
            return;
        end
        
        if sProcess.options.paral.Value  
            poolobj = gcp('nocreate');
            if isempty(poolobj)
                parpool;
            end
        end
           
        path_to_save = workDir;
        DataMat = in_bst_data(DataFile);
        sFile = DataMat.F;

        % Separate the file to length
        isegment = 1;
        nsegment_max = 0;
        
        nget_samples = sProcess.options.binsize.Value{1}*(10^6);   % This is the length of the segment that will be appended.
                                                                   % The reason why we need this is for
                                                                   % machines that don't have enough RAM. The
                                                                   % larger this number, the faster it saves
                                                                   % the files.
        
        while nsegment_max<sFile.prop.samples(2)
            nsegment_min = (isegment-1)*nget_samples;
            nsegment_max = isegment*nget_samples - 1;
            if nsegment_max>sFile.prop.samples(2)
                nsegment_max = sFile.prop.samples(2);
            end
                
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%                     
%             F = in_fread_blackrock(sFile, [nsegment_min, nsegment_max], 1:sFile.header.ChannelCount); % THIS IS JUST FOR BLACKROCK - GENERALIZE
%             F = raw_binary_file.*(10^6); % in_fread_blackrock converts from uV to V. I need it back to uV. 
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 
%             % in_fread_blackrock seems to have the wrong gain. Contact Francois about that. 
%             % signals_and_headers(run).Data.F = signals_and_headers(run).Data.F.*(signals_and_headers(run).extendedheader(1).MaxAnalogValue/signals_and_headers(run).extendedheader(1).MaxDigitalValue); % Convert to uV 

            [F, ~] = in_fread(sFile, ChannelMat, [], [nsegment_min,nsegment_max], [], []);
%           [F, TimeVector] = in_fread(sFile, ChannelMat, iEpoch, SamplesBounds, iChannels, ImportOptions)
            
            
            for ielectrode = 1:sFile.header.ChannelCount
                single_electrode_data = F(ielectrode,:);
                write_to_binary(single_electrode_data, sFile, ielectrode, path_to_save)  % 125 seconds - single core. Waveclus also reads .txt but it took forever to save it to that format
            end        
            isegment = isegment + 1;
            clear raw_binary_file
        end
        disp('')
        disp('')
        disp(['Time for separating channels: ' num2str(toc)])
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%    Convert binaries 2 mat Files    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        bst_progress('start', 'WaveClus spike-sorting', 'Converting binary files to .mat files...');
        if sProcess.options.paral.Value  
            tic
            % This part converts the bin files to mat that is used by WaveClus
            parfor ielectrode = 1:sFile.header.ChannelCount
                convert2mat(sFile, path_to_save, ielectrode) % Single core:18 minutes - parallel: 3.7 minutes (This for a 10 minutes recording)
            end
            
        else
            tic
            % This part converts the bin files to mat that is used by WaveClus
            for ielectrode = 1:sFile.header.ChannelCount
                convert2mat(sFile, path_to_save, ielectrode) % Single core:18 minutes - parallel: 3.7 minutes (This for a 10 minutes recording)
            end
        end
            
            
        disp('')
        disp('')
        disp(['Time for converting to .mat: ' num2str(toc)]) % 271 seconds



        %%%%%%%%%%%%%%%%%%%%%%%%%%% Start the spike sorting %%%%%%%%%%%%%%%%%%%
        bst_progress('start', 'WaveClus spike-sorting', 'Spike-sorting...');
        previous_directory = pwd;
        % The Get_spikes saves the _spikes files at the current directory.

        tic
        all_filenames_in_the_folder = dir(path_to_save);
        cd(all_filenames_in_the_folder(3).folder)
        
        if sProcess.options.paral.Value  
            parfor ielectrode = 1:sFile.header.ChannelCount
                Get_spikes(all_filenames_in_the_folder(ielectrode+2).name)
            end
        else
            for ielectrode = 1:sFile.header.ChannelCount
                Get_spikes(all_filenames_in_the_folder(ielectrode+2).name)
            end
        end
            
            
        disp('')
        disp('')
        disp(['Time for spike-sorting: ' num2str(toc)]) % 528 sec


        %%%%%%%%%%%%%%%%%%%%%% Do the clustering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        bst_progress('start', 'WaveClus spike-sorting', 'Clustering detected spikes...');

        all_filenames_in_the_folder = dir(path_to_save); % This was updated since Get_spikes creates new files

        indices_of_spikes_files = [];
        for iFile = 1:length(all_filenames_in_the_folder)
            if strfind(all_filenames_in_the_folder(iFile).name, '_spikes')
                indices_of_spikes_files = [indices_of_spikes_files iFile];
            end
        end

        % Parfor can't be used on Do_clustering since it creates a file that might be used by several cores
        %     parfor iFile = 1:length(indices_of_spikes_files)
        %         Do_clustering(all_filenames_in_the_folder(indices_of_spikes_files(iFile)).name)
        %     end

        tic
        
        % The optional inputs in Do_clustering have to be true or false, not 1 or 0
        if sProcess.options.paral.Value
            parallel = true;
        else
            parallel = false;
        end
        if sProcess.options.make_plots.Value
            make_plots = true;
        else
            make_plots = false;
        end
        
        Do_clustering(sFile.header.ChannelID, 'parallel',parallel, 'make_plots', make_plots) % Do the clustering in parallel (448 seconds before creating figures)
        disp('')
        disp('')

        disp(['Time for clustering: ' num2str(toc)]) %408 seconds
        
        
        
        %%%%%%%%%%%%%%%%%%%%%  Create Brainstorm Events %%%%%%%%%%%%%%%%%%%
        bst_progress('start', 'WaveClus spike-sorting', 'Saving events file...');
        convert2BrainstormEvents(sFile)
        
        bst_progress('inc', 1);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%   Prepare to exit    %%%%%%%%%%%%%%%%%%%%%%%
    % Turn off parallel processing and return to the initial directory

    if sProcess.options.paral.Value
        if ~isempty(poolobj)
            delete(poolobj);
        end
    end

    cd(previous_directory)
    
    
end



function write_to_binary(data_single_electrode, sFile, ielectrode, path_to_save)
    fid = fopen([path_to_save '\raw_elec' num2str(sFile.header.ChannelID(ielectrode)) '.bin'], 'a');  % THIS JUST APPENDS - CHECK IF IT EXISTS, BEFORE ENTERING THIS LOOP
    
    
    %% IMPORTANT: 10^6 IS NEEDED FOR BLACKROCK: CONVERT TO INTEGERS BEFORE SAVING
    % So I don't use float that would require more space for each file.
    % Check the int16 as well precision used as well. Some acquisition systems might have more precision.
    % MAYBE CONSIDER USING FLOAT.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    data_single_electrode = data_single_electrode*10^6;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fwrite(fid, data_single_electrode,'int16');
    fclose(fid);
    
end

function convert2mat(sFile, path_to_save, ielectrode) % WaveClus doesn't read bin files, I convert it to .mat
    
    fid = fopen([path_to_save '\raw_elec' num2str(sFile.header.ChannelID(ielectrode)) '.bin'], 'rb');
    data = fread(fid, 'int16');
    fclose(fid);  

    %% IMPORTANT: 10^-6 IS NEEDED FOR BLACKROCK: CONVERT TO uV
    % Check if other systems need a conversion from V to uV as well
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    data = data*10^(-6);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sr = sFile.prop.sfreq;
    
    file_delete([path_to_save '\raw_elec' num2str(sFile.header.ChannelID(ielectrode)) '.bin'], 1, 3);
    save([path_to_save '\raw_elec' num2str(sFile.header.ChannelID(ielectrode)) '.mat'],'data','sr') % Waveclus needs variable names: data and sr. No other names work.
    
end



function convert2BrainstormEvents(sFile)

    events = struct;
    events(2).label = [];
    events(2).epochs = [];
    events(2).times = [];
    events(2).color = [];
    events(2).samples = [];
    events(2).reactTimes = [];
    events(2).select = [];
    index = 0;
    for ielectrode = sFile.header.ChannelID'
        
        try
            load(['times_raw_elec' num2str(ielectrode) '.mat'],'cluster_class') % This will fail if the electrode picked up less than 16 spikes. Consider putting a try-catch block

            nNeurons = unique(cluster_class(cluster_class(:,1)>0,1)); % This gives the number of neurons that are picked up on that electrode
            if length(nNeurons)==1
                index = index+1;

                % Write the packet to events
                events(index).label       = ['Spikes Electrode ' num2str(ielectrode)];
                events(index).color       = [rand(1,1),rand(1,1),rand(1,1)];
                events(index).epochs      = ones(1,sum(cluster_class(:,1)~=0));
                events(index).times       = cluster_class(cluster_class(:,1)~=0,2)'./1000; % The timestamps in the cluster_class are in ms
                events(index).samples     = events(index).times.*sFile.prop.sfreq;
                events(index).reactTimes  = [];
                events(index).select      = 1;

            elseif length(nNeurons)>1
                for ineuron = 1:length(nNeurons)
                    % Write the packet to events
                    index = index+1;
                    events(index).label = ['Spikes Electrode ' num2str(ielectrode) ' |' num2str(ineuron) '|'];

                    events(index).color       = [rand(1,1),rand(1,1),rand(1,1)];
                    events(index).epochs      = ones(1,length(cluster_class(cluster_class(:,1)==ineuron,1)));
                    events(index).times       = cluster_class(cluster_class(:,1)==ineuron,2)'./1000; % The timestamps in the cluster_class are in ms
                    events(index).samples     = events(index).times.*sFile.prop.sfreq;
                    events(index).reactTimes  = [];
                    events(index).select      = 1;
                end
            elseif length(nNeurons) == 0
                disp(['Electrode: ' num2str(ielectrode) ' just picked up noise'])
                continue % This electrode just picked up noise
            end
            
        catch
            disp(['Electrode: ' num2str(ielectrode) ' had no clustered spikes'])
        end
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    save('events_UNSUPERVISED.mat','events')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



end

