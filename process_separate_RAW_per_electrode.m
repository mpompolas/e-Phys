function varargout = process_separate_RAW_per_electrode( varargin )

% @=============================================================================

% This process separates the initial raw signal to nChannels binary signals.
% This is performed since most files are too big to load into memory and be
% spike sorted.
% The input is the number of samples that sequentially will be added to
% each file. The larger the number, the faster the procedure. Beware of
% memory issues for a large number.

% Konstantinos Nasiotis
% @=============================================================================






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
    sProcess.Comment     = 'Separate Raw';
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
    sProcess.options.binsize.Value   = {10, 'million samples', 1};
    
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
    
    % Convert all the files in input
    for i = 1:length(sInputs)
        DataFile = file_fullpath(sInputs(i).FileName); % E:\brainstorm_db\Playground\data\Monkey\@rawtest_LFP_EYE\data_0raw_test_LFP_EYE.mat
                                                                                          
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                            
        % MARTIN        
        
        path_to_save = ['E:/brainstorm_db/Playground/data/Monkey/' sTargetStudy.Name]; % I just want the path so I can save the files:
                                                                                       % E:\brainstorm_db\Playground\data\Monkey\@rawtest_LFP_EYE\
                                                                                       % Compare it to the DataFile right above 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                   
                                                                                
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
                                
                                                       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        % MARTIN        

        % ADD SOMETHING HERE TO CHECK IF THE FILES HAVE ALREADY BEEN
        % CREATED. MAYBE CHECK FOR THE EXISTENCE OF THE LAST FILE:
        % raw_elec1003.bin. Without this check, the files would continue getting
        % appended every time we call the function.
        % IF IT EXISTS JUST SHOW A NOTIFICATION TO THE USER THAT IS ALREADY
        % COMPUTED
        
        if exist([path_to_save '/raw_elec' num2str(sFile.header.ChannelID(sFile.header.ChannelCount)) '.bin'])
            STOP
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        
        
        while nsegment_max<sFile.prop.samples(2)
            nsegment_min = (isegment-1)*nget_samples;
            nsegment_max = isegment*nget_samples - 1;
            if nsegment_max>sFile.prop.samples(2)
                nsegment_max = sFile.prop.samples(2);
            end
                
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            % MARTIN              
            raw_binary_file = in_fread_blackrock(sFile, [nsegment_min, nsegment_max], 1:sFile.header.ChannelCount); % THIS IS JUST FOR BLACKROCK - GENERALIZE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
           
            raw_binary_file = raw_binary_file.*(10^6); % in_fread_blackrock converts from uV to V. I need it back to uV. 
            
            % in_fread_blackrock seems to have the wrong gain. Contact
            % Francois about that. 
            % signals_and_headers(run).Data.F = signals_and_headers(run).Data.F.*(signals_and_headers(run).extendedheader(1).MaxAnalogValue/signals_and_headers(run).extendedheader(1).MaxDigitalValue); % Convert to uV 

            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%    REMEMBER THIS     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                 
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % IF RUN OUT OF MEMORY, DELETE ALL THE CREATED BINARY FILES %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            p = gcp('nocreate');
            if isempty(p)
                parpool;
            end
            parfor ielectrode = 1:sFile.header.ChannelCount
                write_to_binary(raw_binary_file(ielectrode,:), sFile, ielectrode, path_to_save)
            end        
            isegment = isegment + 1;
            clear raw_binary_file
        end
    end
end




function write_to_binary(data_single_electrode, sFile, ielectrode, path_to_save)

    ftemp = fopen([path_to_save '/raw_elec' num2str(sFile.header.ChannelID(ielectrode)) '.bin'], 'a');  % THIS JUST APPENDS - CHECK IF IT EXISTS, BEFORE ENTERING THIS LOOP
    fwrite(ftemp, data_single_electrode,'int16');
    fclose(ftemp);    
    
end





