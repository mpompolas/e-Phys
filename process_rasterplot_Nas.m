function varargout = process_timefreq( varargin )
% PROCESS_TIMEFREQ: Computes the time frequency decomposition of any signal in the database.
% 
% USAGE:  sProcess = process_timefreq('GetDescription')
%           sInput = process_timefreq('Run',     sProcess, sInput)
%           TFmask = process_timefreq('GetEdgeEffectMask', Time, Freqs, tfOptions)

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
% Authors: Francois Tadel, 2010-2016

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'e-Phys Raster Plot';
    sProcess.FileTag     = 'raster';
    sProcess.Category    = 'custom';
    sProcess.SubGroup    = 'Exctract';
    sProcess.Index       = 1505;
    sProcess.Description = 'www.in.gr';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data'};
    sProcess.OutputTypes = {'timefreq', 'timefreq', 'timefreq'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    % Options: Sensor types
    sProcess.options.sensortypes.Comment = 'Sensor types or names (empty=all): ';
    sProcess.options.sensortypes.Type    = 'text';
    sProcess.options.sensortypes.Value   = 'MEG, EEG';
    sProcess.options.sensortypes.InputTypes = {'data'};
    sProcess.options.sensortypes.Group   = 'input';
%     % Options: Scouts
%     sProcess.options.clusters.Comment = '';
%     sProcess.options.clusters.Type    = 'scout_confirm';
%     sProcess.options.clusters.Value   = {};
%     sProcess.options.clusters.InputTypes = {'results'};
%     sProcess.options.clusters.Group   = 'input';
%     % Options: Scout function
%     sProcess.options.scoutfunc.Comment    = {'Mean', 'Max', 'PCA', 'Std', 'All', 'Scout function:'};
%     sProcess.options.scoutfunc.Type       = 'radio_line';
%     sProcess.options.scoutfunc.Value      = 1;
%     sProcess.options.scoutfunc.InputTypes = {'results'};
%     sProcess.options.scoutfunc.Group   = 'input';
%     % Options: Time-freq
%     sProcess.options.edit.Comment = {'panel_timefreq_options', 'Morlet wavelet options: '};
%     sProcess.options.edit.Type    = 'editpref';
%     sProcess.options.edit.Value   = [];
%     % Options: Normalize
%     sProcess.options.labelnorm.Comment = '<BR>Spectral flattening:';
%     sProcess.options.labelnorm.Type    = 'label';
%     sProcess.options.normalize.Comment = {'<B>None</B>: Save non-standardized time-frequency maps', '<B>1/f compensation</B>: Multiply output values by frequency'; ...
%                                           'none', 'multiply'};
%     sProcess.options.normalize.Type    = 'radio_label';
%     sProcess.options.normalize.Value   = 'none';
    % Options: Bin size
    sProcess.options.binsize.Comment = 'Bin size: ';
    sProcess.options.binsize.Type    = 'value';
    sProcess.options.binsize.Value   = {0.05, 'ms', 1};
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
    
    % Bin size
    if isfield(sProcess.options, 'binsize') && ~isempty(sProcess.options.binsize) && ~isempty(sProcess.options.binsize.Value) && iscell(sProcess.options.binsize.Value) && sProcess.options.binsize.Value{1} > 0
        bin_size = sProcess.options.binsize.Value{1};
    else
        bst_report('Error', sProcess, sInputs, 'Positive bin size required.');
        return;
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
    
    % === START COMPUTATION ===
    sampling_rate = abs(1. / (tfOPTIONS.TimeVector(2) - tfOPTIONS.TimeVector(1)));
    
    [temp, ~] = in_bst(sInputs(1).FileName);
    nElectrodes = size(temp.ChannelFlag,1); 
    nTrials = length(sInputs);
    nBins = ceil(length(tfOPTIONS.TimeVector) / (bin_size * sampling_rate));
    raster = zeros(nElectrodes, nTrials, nBins);
    
    sample_radius = abs((temp.Time(1) - temp.Time(2)) / 10); % This is used just to extend the first and last bin outside of the Time limits. 
                                                             % If a spike occurs exaclty at the first or last sample, the code crashed. This takes care of that.
    bins = linspace(temp.Time(1) - sample_radius, temp.Time(end) + sample_radius, nBins);
    
    for ifile = 1:length(sInputs)
        [trial, ~] = in_bst(sInputs(ifile).FileName);
        single_file_binning = zeros(nElectrodes, nBins);
        
        for ielectrode = 1:size(trial.F,1)
            for ievent = 1:size(trial.Events,2)
                if strcmp(trial.Events(ievent).label, ['Spikes Electrode ' num2str(ielectrode)])
                     [~, ~, bin_it_belongs_to] = histcounts(trial.Events(ievent).times, bins);
                     unique_bin = unique(bin_it_belongs_to);
                     occurences = [unique_bin; histc(bin_it_belongs_to, unique_bin)];
                     
                     single_file_binning(ielectrode,occurences(1,:)) = occurences(2,:);
                end
            end
            
        end
        
        raster(:, ifile, :) = single_file_binning;
    end
    
    % Get channel file
    sChannel = bst_get('ChannelForStudy', iStudy);
    % Load channel file
    ChannelMat = in_bst_channel(sChannel.FileName);
    tfOPTIONS.ParentFiles = {sInputs.FileName};
    
    % Prepare output file structure
    FileMat.TF = permute(raster, [1 3 2]);
    FileMat.Time = bins;
    FileMat.TFmask = true(size(raster, 2), size(raster, 3));
    FileMat.Freqs = 1:size(FileMat.TF, 3);
    FileMat.Std = [];
    FileMat.Comment = 'Raster Plot';
    FileMat.DataType = 'data';
    FileMat.TimeBands = [];
    FileMat.RefRowNames = [];
    FileMat.RowNames = {ChannelMat.Channel.Name};
    FileMat.Measure = 'power';
    FileMat.Method = 'morlet';
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
    
    % Get output study
    sTargetStudy = bst_get('Study', iStudy);
    % Output filename
    FileName = bst_process('GetNewFilename', bst_fileparts(sTargetStudy.FileName), 'timefreq_rasterplot');
    OutputFiles = {FileName};
    % Save output file and add to database
    bst_save(FileName, FileMat, 'v6');
    db_add_data(tfOPTIONS.iTargetStudy, FileName, FileMat);
    % Display report to user
    bst_report('Info', sProcess, sInputs, 'Success');
    disp('BST> process_timefreq: Success');
end




