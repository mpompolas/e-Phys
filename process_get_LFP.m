function varargout = process_get_LFP( varargin )

% @=============================================================================

% This process loads the binary files created from
% process_separate_RAW_per_electrode and low-passes them at 0.5-300Hz.
% Consecutively, it downsamples to 1000 Hz and converts all of these files to a .ns2 file.

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
    sProcess.Comment     = 'Get LFP';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Import recordings';
    sProcess.Index       = 1201;
    sProcess.Description = 'www.in.gr';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'raw'};
    sProcess.OutputTypes = {'raw'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.isSeparator = 1;
%     sProcess.options.binsize.Comment = 'This will create an LFP.ns2 file';
%     sProcess.options.binsize.Type    = 'value';
% %     sProcess.options.binsize.Value   = {10, 'million samples', 1};
    
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
    
    % Get channel file
    sChannel = bst_get('ChannelForStudy', iStudy);
    % Load channel file
    ChannelMat = in_bst_channel(sChannel.FileName);
    
    
    
    % Convert all the files in input
    for i = 1:length(sInputs)
        DataFile = file_fullpath(sInputs(i).FileName); % E:\brainstorm_db\Playground\data\Monkey\@rawtest_LFP_EYE\data_0raw_test_LFP_EYE.mat
                                                                                          
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                            
        % MARTIN        
        
        path_to_load = ['E:/brainstorm_db/Playground/data/Monkey/' sTargetStudy.Name];   % I just want the path so I can load the files:
                                                                                         % E:\brainstorm_db\Playground\data\Monkey\@rawtest_LFP_EYE\
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                   
                                                                                
        DataMat = in_bst_data(DataFile);
        sFile = DataMat.F;
        
        % Add the header to the file
        add_header(sFile, path_to_load, ChannelMat)
        
        % Add the LFP to the file
        add_lfp(sFile, path_to_load)
        
    end
end







function add_header(sFile, path_to_load, ChannelMat)
    % Copy the header from the initial file and then make changes on it
    % where needed.
    
    finitial = fopen(sFile.filename,'rb');
    initial_header_extended_data = fread(finitial,314 + 66 * sFile.header.ChannelCount + 9); % All the headers together from that file: header, extended header, data header.
    fclose(finitial);

    % Copy that header to the new LFP file and change the Labels of the
    % electrodes to "LFP ichannel"
    
    fid = fopen([path_to_load '/' sFile.comment(1:end-4) '_LFP.ns2'],'w'); %*.ns2:  Continuous LFP data sampled at 1 KHz
    fwrite(fid,initial_header_extended_data);
    
    fseek(fid, 14, 'bof');
    fwrite(fid,'1 ksamp/sec     ','uint8'); % Change the sampling rate on the header.
    
    fseek(fid, 290, 'bof');
    fwrite(fid,1000,'uint32'); % Change the sampling rate on the header.
    
    fseek(fid, 314 + 66*(sFile.header.ChannelCount) + 5, 'bof');
    fwrite(fid, ceil(sFile.header.DataPoints/(sFile.header.SamplingFreq/1000)),'uint32'); % Change the number of samples to the about to be filtered signal's length
    
    for ielectrode = 1:sFile.header.ChannelCount
        if strfind(ChannelMat.Channel(ielectrode).Name,'raw')
            fseek(fid,314 + 4 + 66*(ielectrode-1),'bof');
            fwrite(fid,'LFP');
        end
    end
    fclose(fid);
end
   



function add_lfp(sFile, path_to_load)

    LFP = zeros(sFile.header.ChannelCount, ceil(sFile.header.DataPoints/(sFile.header.SamplingFreq/1000)),'int16');
    
    % These are small files, it can be done in parallel.
    p = gcp('nocreate');
    if isempty(p)
        parpool;
    end
    parfor ielectrode = 1:sFile.header.ChannelCount
        data_filtered = filter_and_downsample_files(sFile, path_to_load, ielectrode)
        LFP(ielectrode,:) = data_filtered;
    end  
    
    % Convert back to .ns2
    % *.ns2:  Continuous LFP data sampled at 1 KHz
    
    ffinal = fopen([path_to_load '/' sFile.comment(1:end-4) '_LFP.ns2'],'a');

    fwrite(ffinal, LFP,'int16');
    fclose(ffinal);
    
    path_to_load = strrep(path_to_load,'/','\'); % The problem was using '/' instead of '\'
    
    OutputFiles = import_raw({[path_to_load '\' sFile.comment(1:end-4) '_LFP.ns2']}, 'EEG-RIPPLE', 2); 
%     OutputFiles = import_raw({[path_to_load '/LFP.ns2']}, 'EEG-RIPPLE', iSubject, ImportOptions)

%-------------------------------------
% Martin 
% Get the iSubject value to add above.
% ImportOptions can be ommitted.
%-------------------------------------
    
    
  
end




function data_filtered = filter_and_downsample_files(sFile, path_to_load, ielectrode)
    disp(num2str(ielectrode))

    fid = fopen([path_to_load '/raw_elec' num2str(sFile.header.ChannelID(ielectrode)) '.bin'], 'r');    
    data = fread(fid, 'int16')';
    fclose(fid);

    [data_filtered, FiltSpec, Messages] = bst_bandpass_hfilter(data, sFile.header.SamplingFreq, 0.5, 300, 0, 0);
%         [data_filtered, FiltSpec, Messages] = bst_bandpass_hfilter(data, sFile.header.SamplingFreq, HighPass, LowPass, isMirror, isRelax);

    data_filtered = downsample(data_filtered,sFile.header.SamplingFreq/1000);  % The file now has a different sampling rate (fs/30) = 1000Hz. This has to be stored somewhere

end






