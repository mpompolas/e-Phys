function varargout = process_spike_field_coherence( varargin )
% PROCESS_SPIKE_FIELD_COHERENCE: Computes the spike field coherence.
% 
% USAGE:    sProcess = process_spike_field_coherence('GetDescription')
%        OutputFiles = process_spike_field_coherence('Run', sProcess, sInput)

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
% Authors: Konstantinos Nasiotis, 2017; Martin Cousineau, 2017

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Spike Field Coherence';
    sProcess.FileTag     = 'SFC';
    sProcess.Category    = 'custom';
    sProcess.SubGroup    = 'Exctract';
    sProcess.Index       = 2506;
    sProcess.Description = 'www.in.gr';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data'};
    sProcess.OutputTypes = {'timefreq'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    % Options: Sensor types
    sProcess.options.sensortypes.Comment = 'Sensor types or names (empty=all): ';
    sProcess.options.sensortypes.Type    = 'text';
    sProcess.options.sensortypes.Value   = 'EEG';
    sProcess.options.sensortypes.InputTypes = {'data'};
    sProcess.options.sensortypes.Group   = 'input';
    % Options: Bin size
    sProcess.options.timewindow.Comment  = 'Spike Time window     :';
    sProcess.options.timewindow.Type     = 'timewindow';
    sProcess.options.timewindow.Value    = [];
   
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    % Initialize returned values
    OutputFiles = {};
    % Extract method name from the process name
    strProcess = strrep(strrep(func2str(sProcess.Function), 'process_', ''), 'timefreq', 'morlet');
    
    % Add other options
    tfOPTIONS.Method = strProcess;
    if isfield(sProcess.options, 'sensortypes')
        tfOPTIONS.SensorTypes = sProcess.options.sensortypes.Value;
    else
        tfOPTIONS.SensorTypes = [];
    end    
    
    % If a time window was specified
    if isfield(sProcess.options, 'timewindow') && ~isempty(sProcess.options.timewindow) && ~isempty(sProcess.options.timewindow.Value) && iscell(sProcess.options.timewindow.Value)
        tfOPTIONS.TimeWindow = sProcess.options.timewindow.Value{1};
    elseif ~isfield(tfOPTIONS, 'TimeWindow')
        tfOPTIONS.TimeWindow = [];
    end
    
    % Output
    if isfield(sProcess.options, 'avgoutput') && ~isempty(sProcess.options.avgoutput) && ~isempty(sProcess.options.avgoutput.Value)
        if sProcess.options.avgoutput.Value
            tfOPTIONS.Output = 'average';
        else
            tfOPTIONS.Output = 'all';
        end
    end
    
    tfOPTIONS.TimeVector = in_bst(sInputs(1).FileName, 'Time');

    
    % === OUTPUT STUDY ===
    % Get output study
    [~, iStudy, ~] = bst_process('GetOutputStudy', sProcess, sInputs);
    tfOPTIONS.iTargetStudy = iStudy;
   
    % Get channel file
    sChannel = bst_get('ChannelForStudy', iStudy);
    % Load channel file
    ChannelMat = in_bst_channel(sChannel.FileName);
    
    
    % === START COMPUTATION ===
    sampling_rate = round(abs(1. / (tfOPTIONS.TimeVector(2) - tfOPTIONS.TimeVector(1))));
    
    [temp, ~] = in_bst(sInputs(1).FileName);
    
    nElectrodes = size(temp.ChannelFlag,1); 
    nTrials = length(sInputs);
    
    all_trials_FFTs = cell(nElectrodes, 2);
    
    % Optimize this
    for ifile = 1:nTrials
        [trial, ~] = in_bst(sInputs(ifile).FileName);
        
        for ielectrode = 1:nElectrodes
            disp(ielectrode)
            for ievent = 1:size(trial.Events, 2)
                if strcmp(trial.Events(ievent).label, ['Spikes Electrode ' num2str(ielectrode)])
                    
                    events_within_segment = trial.Events(ievent).samples(trial.Events(ievent).times > trial.Time(1) + abs(sProcess.options.timewindow.Value{1}(1)) & ...
                                                                         trial.Events(ievent).times < trial.Time(end) - abs(sProcess.options.timewindow.Value{1}(2)));
                    
                    for ispike = 1:length(events_within_segment)
                        all_trials_FFTs{ielectrode,1} = [all_trials_FFTs{ielectrode,1};
                                                         trial.F(ielectrode, ...
                                                                 round(length(trial.Time) / 2) + events_within_segment(ispike) - abs(sProcess.options.timewindow.Value{1}(1)) * sampling_rate: ...
                                                                 round(length(trial.Time) / 2) + events_within_segment(ispike) + abs(sProcess.options.timewindow.Value{1}(2)) * sampling_rate ...
                                                                ) ...
                                                        ];
                                                                  
                        [temp_FFT, Freqs] = compute_FFT(all_trials_FFTs{ielectrode,1}(end,:), trial.Time);  
                        all_trials_FFTs{ielectrode,2} = [all_trials_FFTs{ielectrode,2} ; temp_FFT];
                    end
                   
                end
            end
            
        end
    end 
        
    
    
    % Calculate the SFC
    SFC = zeros(nElectrodes, 1, ceil(length(temp.Time) / 2)); %MAKE SURE THE SECOND DIMENSION IS CORRECT HERE

    
    
    for ielectrode = 1:nElectrodes
        
        single_electrode_FFTs   = abs(all_trials_FFTs{ielectrode,2}) .^ 2; % Power
        single_electrode_trials = all_trials_FFTs{ielectrode,1};
        

        [FFT_of_average, ~]   = compute_FFT(sum(single_electrode_trials), trial.Time);  
        SFC(ielectrode, 1, :) = (abs(FFT_of_average) .^ 2) ./ sum(single_electrode_FFTs);
    
    end
    
    SFC(isnan(SFC)) = 0; % Convert NaN to zeros cause Brainstorm won't read it otherwises
    
    
    
    tfOPTIONS.ParentFiles = {sInputs.FileName};
    
    % Prepare output file structure
    FileMat.TF = SFC;
    FileMat.Time = [trial.Time(1), trial.Time(end)];
    FileMat.TFmask = [];
    FileMat.Freqs = Freqs;
    FileMat.Std = [];
    FileMat.Comment = 'Spike Field Coherence Plot';
    FileMat.DataType = 'data';
    FileMat.TimeBands = [];
    FileMat.RefRowNames = [];
    FileMat.RowNames = {ChannelMat.Channel.Name};
    FileMat.Measure = 'power';
    FileMat.Method = 'fft';
    FileMat.DataFile = []; % Leave blank because multiple parents
    FileMat.SurfaceFile = [];
    FileMat.GridLoc = [];
    FileMat.GridAtlas = [];
    FileMat.Atlas = [];
    FileMat.HeadModelFile = [];
    FileMat.HeadModelType = [];
    FileMat.nAvg = [];
    FileMat.ColormapType = [];
    FileMat.DisplayUnits = [];
    FileMat.Options = tfOPTIONS;
    FileMat.History = [];
    
    
%     % Options
%     FileMat.Options.Comment = 'FFT';
%     FileMat.Options.Method = 'fft'; % Maybe this would work.
%     FileMat.Options.Freqs = [];
%     FileMat.Options.TimeVector = trial.Time;
%     FileMat.Options.TimeBands = [];
%     FileMat.Options.TimeWindow = [trial.Time(1) , trial.Time(end)];
%     FileMat.Options.ClusterFuncTime = 'none';
%     FileMat.Options.Measure = 'power';
%     FileMat.Options.Output = 'all';
%     FileMat.Options.RemoveEvoked = 0;
%     FileMat.Options.MorletFc = 1;
%     FileMat.Options.MorletFwhmTc = 3;
%     FileMat.Options.WinLength = [];
%     FileMat.Options.WinOverlap = 50;
%     FileMat.Options.isMirror = 0;
%     FileMat.Options.SensorTypes = '';
%     FileMat.Options.Clusters = [];
%     FileMat.Options.ScoutFunc = [];
%     FileMat.Options.SurfaceFile = [];
%     FileMat.Options.iTargetStudy = [];
%     FileMat.Options.SaveKernel = 0;
%     FileMat.Options.nComponenets = 1;
%     FileMat.Options.NormalizeFunc = 'none';
    
    
    % Get output study
    sTargetStudy = bst_get('Study', iStudy);
    % Output filename
    FileName = bst_process('GetNewFilename', bst_fileparts(sTargetStudy.FileName), 'timefreq_fft_spike_field_coherence');
    OutputFiles = {FileName};
    % Save output file and add to database
    bst_save(FileName, FileMat, 'v6');
    db_add_data(tfOPTIONS.iTargetStudy, FileName, FileMat);
    % Display report to user
    bst_report('Info', sProcess, sInputs, 'Success');
    disp('BST> process_spike_field_coherence: Success');
end





function [TF, Freqs] = compute_FFT(F, time)

    % Next power of 2 from length of signal
    nTime = length(time);
    % NFFT = 2^nextpow2(nTime);    % Function fft() pads the signal with zeros before computing the FT
    NFFT = nTime;                  % No zero-padding: Nfft = Ntime
    sfreq = 1 / (time(2) - time(1));
    % Positive frequency bins spanned by FFT
    Freqs = sfreq / 2 * linspace(0, 1, NFFT / 2 + 1);
    % Keep only first and last time instants
    time = time([1, end]);
    % Remove mean of the signal
    F = bst_bsxfun(@minus, F, mean(F,2));
    % Apply a hamming window to signal
    F = bst_bsxfun(@times, F, bst_window('hamming', size(F,2))');
    % Compute FFT
    Ffft = fft(F, NFFT, 2);
    % Keep only first half
    % (x2 to recover full power from negative frequencies)
    TF = 2 * Ffft(:, 1:NFFT / 2 + 1) ./ nTime;
    % Permute dimensions: time and frequency
    TF = permute(TF, [1 3 2]);
    
end




