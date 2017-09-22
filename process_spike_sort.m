function varargout = process_spike_sort( varargin )

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
    sProcess.Comment     = 'MEA AUTO Spike Sorter';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Exctract';
    sProcess.Index       = 1200;
    sProcess.Description = 'www.in.gr';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'raw'};
    sProcess.OutputTypes = {'raw'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.isSeparator = 1;
%     sProcess.options.binsize.Comment = 'Samplesize for appending: ';
%     sProcess.options.binsize.Type    = 'value';
%     sProcess.options.binsize.Value   = {10, 'million samples', 1};
    
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
        [path_to_load, filebase, fileext] = bst_fileparts(file_fullpath(sInputs(i).FileName));

        start_spike_sorter(path_to_load, filebase)
        
        
    end
end




function start_spike_sorter(path_to_load, filebase)

        % function batchsort2(par)
    disp('Starting batch sort....');
    load('parameters.mat') % This has to be in the path. I currently just copied it at: C:\Users\McGill\Desktop\Spike sorter

    %Fully automated sorting of extracellular data, optimized for parallel
    %processing, modular code.
    par.pathname = path_to_load;
    par.filename = filebase; % the link filename

    pathname = par.pathname;
    cd (pathname)

    folderFiles = dir;
    iheaders = find(strcmp({folderFiles.name},'channel.mat'));
    electrode_files = nan(length(folderFiles),1);
    electrode_indices = strrep({folderFiles.name},'raw_elec','');
    electrode_indices = strrep(electrode_indices,'.bin','')';

    for ifile = 1:length(folderFiles)
        if strfind(folderFiles(ifile).name,'raw_elec')
            electrode_files(ifile) = str2num(electrode_indices{ifile});
        end
    end

    load(folderFiles(iheaders).name)
    link_file = load([pathname '/' par.filename]);
    clear iheaders

    par.sr=link_file.F.prop.sfreq; % Load the sampling rate of the file
    clear data_info

    handles.par=par;

    handles.file_name=par.filename;%random non-GUI fix
    handles.datatype='ASCII';
    filename_noext=par.filename(7:end-4);%filename without 'data_0' and '.mat'. Brainsotrm saves it with data_0 as a prefix
    filename=par.filename;

    numchan = length(Channel);

    %Find the number of samples in a channel
    chanEl = link_file.F.prop.samples(2)+1;

    %Never process more than 20 minutes at a time - This doesnt make any sense

    % The following splits the code
    numsplits = 1;
    startsplit = 1;
    maxlen = 60*par.sr*20;   %maxlen = 20 minutes = 60*30000*20


    for mm = startsplit:numsplits
        par = handles.par;

        offset = 0;

        %rereference
        refchan=[];
        multipliers = ones(numchan,1);


        % Setup for template matching to previous session
        if strcmp(par.clusterType,'tempmatch')
            %load up data from prev. results. We need the features (nspk) and
            %classifications of the old spikes
            prevResults=load([par.templatePathname 'sorted_' par.templateFilename],'spikes','cluster_class');
            prevSpikes = prevResults.spikes;
            prevClass = prevResults.cluster_class;
            clear prevResults
        end


        channels = [];
        for ichannel = 1:length(Channel)
            index_channel = str2num(strrep(Channel(ichannel).Name,'raw ','')); 
            channels = [channels index_channel];
        end

        %% Spike sort per channel
        for i=1:numchan
            chani=channels(i);

            ichannel = find(electrode_files == chani);

            fid = fopen(folderFiles(ichannel).name,'rb');
            chandata = fread(fid, [1, chanEl], 'int16')';
            fclose(fid);

            Bandpass_low  = 300;  % original  500; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Bandpass_high = 8000; % original 8000;

            %denoise via ref. channel
            Wn = [Bandpass_low,Bandpass_high]/par.sr;
            [ha,hb] = ellip(2,.1,40,Wn);
            chandata = filtfilt(ha,hb,chandata);

            upSampData = chandata;
            clear chandata;

            %Filter data
            disp('Filtering');    


            %Super optimal - This is for filtering (experimental)
            if isfield(handles.par,'optimalf') && handles.par.optimalf
                athresh = handles.par.stdmin;

                offsets = [0,2,2];
                a = load('~/SharedDocuments/Patrick/optimalfilteringpaper/multifilters4.mat');
                wfilters = a.filters;
                indices = [];

                thresh = -norminv((1-normcdf(athresh))/size(wfilters,1));
                handles.par.stdmin = thresh;

                for ii = 1:size(wfilters,1)
                    upSampData2 = conv(upSampData,-wfilters(ii,end:-1:1),'same');
                    [~,thr,index] = amp_detect(upSampData2-median(upSampData2),handles);
                    index = index + offsets(ii);
                    indices = [indices,index];
                    nums(ii) = length(index);
                end

                nums

                handles.par.stdmin = athresh;

                indices = sort(indices);
                badidx = indices < [0,indices(1:end-1)] + 6;
                indices = indices(~badidx);
                badidx = indices < [0,indices(1:end-1)] + 6;
                indices = indices(~badidx);
                badidx = indices < [0,indices(1:end-1)] + 6;
                indices = indices(~badidx);

                spikes = upSampData2(bsxfun(@plus,indices',-7:7));
                [~,optpos] = min(spikes,[],2);
                index = optpos+indices(:)-7;

                badidx = index < [0;index(1:end-1)] + 6;
                index = index(~badidx);

                spikes = upSampData2(bsxfun(@plus,index,-handles.par.w_pre:handles.par.w_post+1));
                upSampData = upSampData2;
            else % This is for filtering (Regular)
                Wn = [Bandpass_low,Bandpass_high]/par.sr;
                [ha,hb] = ellip(2,.1,40,Wn);
                upSampData = filtfilt(ha,hb,upSampData);
                upSampData = upSampData'; % KONSTANTINOS
                upSampData = bpfiltout(upSampData,[3998 4002 4600 4750 7998 8002],par.sr,2);

                sdat = upSampData(1:100:end);
                md = median(sdat(:));
    % % % % % %             handles.par.stdmin = 6; % Added by Theo for Praveen data
                [spikes,thr,index]=amp_detect_nsignals(upSampData-md,upSampData'-md,handles); % spikes holds the spike waveform, index probably when it occured- check this out


            end

            clear sdat sdat2

            thrmax=10*thr;

            fprintf('%d spikes detected\n',size(spikes,1));

            %No use in remembering all this data
            upSampDataSnip = upSampData(1:min(length(upSampData),5e6));
            lenData = length(upSampData)/handles.par.sr;

            index=index-1; %to account for the first sample being at time 0;
            index = index*1e3/handles.par.sr;                     %spike times in ms.

            switch par.clusterType
                case 'spc'

                    clear upSampData;

                    % LOAD SPIKES
                    handles.par.fname = ['data_' filename];   %filename for interaction with SPC
                    nspk = size(spikes,1);
                    naux = min(handles.par.max_spk,size(spikes,1));
                    handles.par.min_clus = max(handles.par.min_clus_abs,handles.par.min_clus_rel*naux);

                    %Patrick: new feature: minimum number of spikes/s
                    if ~isfield(handles.par,'min_clus_rate')
                        handles.par.min_clus_rate = .5;
                    end

                    handles.par.min_clus = max(handles.par.min_clus,round(handles.par.min_clus_rate*lenData));

                    [inspk] = wave_features(spikes,handles);                %Extract spike features.


                    handles.par.fname_in = 'tmp_data';
                    fname_in = handles.par.fname_in;
                    %    save([fname_in],'inspk_aux','-ascii');                         %Input file for SPC
                    handles.par.fnamesave = [handles.par.fname '_' ...
                        filename(1:end-4)];                                %filename if "save clusters" button is pressed
                    handles.par.fname = [handles.par.fname '_wc'];             %Output filename of SPC
                    handles.par.fnamespc = handles.par.fname;

                    % SELECTION OF SPIKES FOR SPC
                    if handles.par.permut == 'n'
                        % GOES FOR TEMPLATE MATCHING IF TOO MANY SPIKES.
                        if size(spikes,1)> handles.par.max_spk;
                            % take first 'handles.par.max_spk' spikes as an input for SPC
                            inspk_aux = inspk(1:naux,:);
                        else
                            inspk_aux = inspk;
                        end

                        disp(['number of spikes' num2str(numelindex)]);
                        %INTERACTION WITH SPC
                        save(handles.par.fname_in,'inspk_aux','-ascii');
                        if numel(index)>11 %minimum number of spikes for spc to work
                            [clu, tree] = run_cluster(handles);
                            [temp] = find_temp(tree,handles);
                        else
                            clu=zeros(21,numel(index)+2);
                            clu(:,1)=0:20;
                            clu(:,2)=0:0.01:0.2;
                            tree=clu;
                            temp=1;

                        end

                        %DEFINE CLUSTERS
                        class1=find(clu(temp,3:end)==0);
                        class2=find(clu(temp,3:end)==1);
                        class3=find(clu(temp,3:end)==2);
                        class4=find(clu(temp,3:end)==3);
                        class5=find(clu(temp,3:end)==4);
                        class0=setdiff(1:size(spikes,1), sort([class1 class2 class3 class4 class5]));
                        whos class*

                    else
                        % GOES FOR TEMPLATE MATCHING IF TOO MANY SPIKES.
                        if size(spikes,1)> handles.par.max_spk;
                            % random selection of spikes for SPC
                            ipermut = randperm(length(inspk));
                            ipermut(naux+1:end) = [];
                            inspk_aux = inspk(ipermut,:);
                        else
                            ipermut = randperm(size(inspk,1));
                            inspk_aux = inspk(ipermut,:);
                        end

                        save(handles.par.fname_in,'inspk_aux','-ascii');
                        if numel(index)>11 %minimum number of spikes for spc binary to work
                            [clu, tree] = run_cluster(handles);
                            [temp] = find_temp(tree,handles);
                        else
                            clu=zeros(21,numel(index)+2);
                            clu(:,1)=0:20;
                            clu(:,2)=0:0.01:0.2;
                            tree=clu;
                            temp=1;

                        end

                        %DEFINE CLUSTERS
                        class1=ipermut(find(clu(temp,3:end)==0));
                        class2=ipermut(find(clu(temp,3:end)==1));
                        class3=ipermut(find(clu(temp,3:end)==2));
                        class4=ipermut(find(clu(temp,3:end)==3));
                        class5=ipermut(find(clu(temp,3:end)==4));
                        class0=setdiff(1:size(spikes,1), sort([class1 class2 class3 class4 class5]));
                        whos class*
                    end



                    % IF TEMPLATE MATCHING WAS DONE, THEN FORCE
                    if (size(spikes,1)> handles.par.max_spk | ...
                            (handles.par.force_auto == 'y'));
                        classes = zeros(size(spikes,1),1);
                        if length(class1)>=handles.par.min_clus; classes(class1) = 1; end
                        if length(class2)>=handles.par.min_clus; classes(class2) = 2; end
                        if length(class3)>=handles.par.min_clus; classes(class3) = 3; end
                        if length(class4)>=handles.par.min_clus; classes(class4) = 4; end
                        if length(class5)>=handles.par.min_clus; classes(class5) = 5; end
                        f_in  = spikes(classes~=0,:);
                        f_out = spikes(classes==0,:);
                        class_in = classes(find(classes~=0),:);
                        class_out = force_membership_wc(f_in, class_in, f_out, handles);
                        classes(classes==0) = class_out;
                        class0=find(classes==0);
                        class1=find(classes==1);
                        class2=find(classes==2);
                        class3=find(classes==3);
                        class4=find(classes==4);
                        class5=find(classes==5);
                    end

                case 'tempmatch'

                    naux = min(handles.par.max_spk,size(spikes,1));
                    handles.par.min_clus = max(handles.par.min_clus_abs,handles.par.min_clus_rel*naux);

                    %Patrick: new feature: minimum number of spikes/s
                    if ~isfield(handles.par,'min_clus_rate')
                        handles.par.min_clus_rate = .5;
                    end

                    handles.par.min_clus = max(handles.par.min_clus,round(handles.par.min_clus_rate*lenData));


                    oldSpikes=prevSpikes{chani};
                    oldCluster_class=prevClass{chani};
                    oldClasses=oldCluster_class(:,1);
                    if size(oldSpikes,1) ~= length(oldClasses) || isempty(oldSpikes)
                        %Edit Patrick: handle some weird race condition
                        %Garbage
                        f_in = zeros(0,size(spikes,2));
                        class_in = zeros(0,size(spikes,2));
                    else
                        f_in=oldSpikes(oldClasses~=0,:);
                        class_in = oldClasses(oldClasses~=0,:);
                    end


                    if max(oldClasses) > 0
                        [~,sortidx] = sort(oldSpikes(:,56));
                        mspk = mean(oldSpikes(sortidx( ... 
                                    round(length(sortidx)*.1):round(length(sortidx)*.9)),:));

                        load wfilter

                        B = mspk*W*W;
                        filteredDat = conv(upSampData,-B(end:-1:1),'same');
                        filteredDat = filteredDat - median(filteredDat);
                        filteredDat = [filteredDat(9:end),zeros(1,8)];

                        [spikes,thr,index1]=amp_detect_nsignals(upSampData,filteredDat'-median(filteredDat(:)),handles);

                        index=index1-1; %to account for the first sample being at time 0;
                        index = index*1e3/handles.par.sr;                     %spike times in ms.
                    end


                    handles.par.template_type = 'nb';
                    f_out=spikes;
                    class_out = force_membership_wc(oldSpikes, oldClasses, f_out, handles);
                    classes = class_out;
                    class0=find(classes==0);
                    class1=find(classes==1);
                    class2=find(classes==2);
                    class3=find(classes==3);
                    class4=find(classes==4);
                    class5=find(classes==5);


                    %create an empty clu/tree/temp
                    clu=zeros(21,numel(index)+2);
                    clu(:,1)=0:20;
                    clu(:,2)=0:0.01:0.2;
                    tree=clu;
                    temp=1;
                    nspk = size(spikes,1);
                    [inspk] = wave_features(spikes,handles);                %Extract spike features.
                    ipermut = randperm(length(inspk));
            end



            clus_pop = [];
            ylimit = [];


            cluster=zeros(nspk,2);
            if ~isempty(index)
                cluster(:,2)= index';
            end

            num_clusters = length(find([length(class1) length(class2) length(class3)...
                length(class4) length(class5) length(class0)] >= handles.par.min_clus));
            clus_pop = [clus_pop length(class0)];

            %Assign things to proper clusters
            for ii = 1:5
                tgt = eval(['class' num2str(ii)]);
                if length(tgt) > handles.par.min_clus
                    cluster(tgt(:),1)=ii;
                end
            end


            %check for erronious clustering (this is a hack)
            a=cluster(:,1)>10;
            cluster(a,1)=0;

            %store results-------------------------------------------------------
            %---------------------------------------------------------------------

            cluster(:,2) = cluster(:,2)+(mm-1)*maxlen/30;


            %correct for upsampling
            handles.par.sr=handles.par.sr;%/handles.par.int_factor;
            par_out{chani}=handles.par;
            cluster_class_out{chani}=cluster;
            tree_out{chani}=tree;
            clu_out{chani}=clu;
            spikes_out{chani}=spikes;
            thr_out{chani}=thr;
            if  handles.par.permut == 'y'
                ipermut_out{chani}=ipermut;
            end
            inspk_out{chani}=inspk;
            temp_out{chani}=temp;

            Signals(chani).Spikes = spikes_out{chani};
            Signals(chani).Cluster_class = cluster_class_out{chani};
            clear spikes cluster_class tree thr inspk clu temp

            disp(['Electrode ' num2str(chani) ' complete!']);



        end


        spikes=spikes_out;
        par=par_out;
        cluster_class=cluster_class_out;
        clu=clu_out;
        tree1=tree_out;
        thr=thr_out;
        inspk=inspk_out;
        temp=temp_out;

        fnamesuffix = '';
        fname = filename;
        if mm > 1
            fname = filename(1:end-4);
            fnamesuffix = sprintf('_%02d.mat',mm);
        end

        if  handles.par.permut == 'y'
            ipermut=ipermut_out;
            %First save locally
            save(['sorted_' fname fnamesuffix], '-v7.3','spikes', 'par', 'cluster_class', ...
                'clu', 'tree1', 'inspk', 'thr', 'ipermut','refchan', 'temp','multipliers');
            save(['original_batch_sorted_' fname fnamesuffix], '-v7.3','spikes', 'par', 'cluster_class', ...
                'clu', 'tree1', 'inspk', 'thr', 'ipermut','refchan', 'temp','multipliers');
        else

            save(['sorted_' fname fnamesuffix], '-v7.3','spikes', 'par', 'cluster_class', ...
                'clu', 'tree1', 'inspk', 'thr', 'refchan', 'temp','multipliers');
            save(['original_batch_sorted_' fname fnamesuffix], '-v7.3','spikes', 'par', 'cluster_class', ...
                'clu', 'tree1', 'inspk', 'thr', 'refchan', 'temp','multipliers');

        end

        %%%%% Save spikes in Brainstorm events format %%%%%%%

        % Check if the events file exists. If not, create one.
        % This means there are no experiment events at all, not just spike-events
        if ~(exist([pathname '/events.mat'], 'file') == 2)
            events = struct;
            events(2).label = [];
            events(2).epochs = [];
            events(2).times = [];
            events(2).color = [];
            events(2).samples = [];
            events(2).reactTimes = [];
            events(2).select = [];
            save([pathname '/events.mat'],'events')
            clear events
        end


        load([pathname '/events.mat']) % These are the events from the acquisition system
        ientries_to_delete = [];
        for ielectrode = channels

            index = find(strcmp({events.label},['Spikes Electrode ' num2str(ielectrode)]));

            if ~isempty(cluster_class(1,ielectrode)) % The number of the electrode might not be the same as the socket on the acquisition system, so it would be empty
                if isempty(index)
                    index = size(events,2)+1;
                    events(index).label = ['Spikes Electrode ' num2str(ielectrode)];
                end

                nNeurons = unique(cluster_class{1,ielectrode}(cluster_class{1,ielectrode}(:,1)>0,1)); % This gives the number of neurons that are picked up on that electrode
                if length(nNeurons)==1

                    % Write the packet to events

                    events(index).color           = [rand(1,1),rand(1,1),rand(1,1)];
                    events(index).epochs          = ones(1,sum(cluster_class{1,ielectrode}(:,1)~=0));
                    events(index).times           = cluster_class{1,ielectrode}(cluster_class{1,ielectrode}(:,1)~=0,2)'./1000; % The timestamps in the cluster_class are in ms
                    events(index).samples         = events(index).times.*par{1,1}.sr;
                    events(index).reactTimes      = [];
                    events(index).select          = 1;

                elseif length(nNeurons)>1
                    ientries_to_delete = [ientries_to_delete index];

                    for ineuron = 1:length(nNeurons)
                        % Write the packet to events
                        index = size(events,2)+1;
                        events(index).label = ['Spikes Electrode ' num2str(ielectrode) ' |' num2str(ineuron) '|'];

                        events(index).color       = [rand(1,1),rand(1,1),rand(1,1)];
                        events(index).epochs      = ones(1,length(cluster_class{1,ielectrode}(cluster_class{1,ielectrode}(:,1)==ineuron,1)));
                        events(index).times       = cluster_class{1,ielectrode}(cluster_class{1,ielectrode}(:,1)==ineuron,2)'./1000; % The timestamps in the cluster_class are in ms
                        events(index).samples     = events(index).times.*par{1,1}.sr;
                        events(index).reactTimes  = [];
                        events(index).select      = 1;
                    end

                elseif length(nNeurons)==0   

                    events(index).epochs          = [];
                    events(index).times           = [];
                    events(index).samples         = [];

                    ientries_to_delete = [ientries_to_delete index];

                end
            end

        end
        events(ientries_to_delete) = []; % Delete the entry. This will be substituted by 2-3 new entries, based on the number of neurons that are picked up from that electrode.

        % Get rid of the empty entries and order the events alphabetically based on
        % their label.
        select_events_that_have_values = false(1,length(events));
        for i = 1:length(events)
            select_events_that_have_values(i) = ~isempty(events(i).label);
        end
        events = events(select_events_that_have_values);

        [~, order_alphabetically_indices] = sort(string({events.label}));
        events = events(order_alphabetically_indices);
        if isfield(events,'UnitClassificationNumber') && isfield(events,'indicesOnPacketsStruct')
            events = rmfield(events,{'UnitClassificationNumber';'indicesOnPacketsStruct'});
        elseif isfield(events,'UnitClassificationNumber') && isfield(events,'PacketInsertionReason')
            events = rmfield(events,{'UnitClassificationNumber';'PacketInsertionReason'});
        end

        save([pathname '/events_original_batch_sorted_',  filename_noext ,'.mat'],'events', 'ipermut', 'spikes', 'par', 'cluster_class', 'clu', 'tree1', 'inspk', 'thr', 'refchan', 'temp','multipliers');
        save([pathname '/events_BRAINSTORM_AUTO_sorted_', filename_noext ,'.mat'],'events')
        disp('Unsupervised sorting is complete.');

end


    
end





