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
    sProcess.InputTypes  = {'data', 'results', 'matrix'};
    sProcess.OutputTypes = {'timefreq', 'timefreq', 'timefreq'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    % Options: Sensor types
    sProcess.options.sensortypes.Comment = 'Sensor types or names (empty=all): ';
    sProcess.options.sensortypes.Type    = 'text';
    sProcess.options.sensortypes.Value   = 'MEG, EEG';
    sProcess.options.sensortypes.InputTypes = {'data'};
    sProcess.options.sensortypes.Group   = 'input';
    % Options: Scouts
    sProcess.options.clusters.Comment = '';
    sProcess.options.clusters.Type    = 'scout_confirm';
    sProcess.options.clusters.Value   = {};
    sProcess.options.clusters.InputTypes = {'results'};
    sProcess.options.clusters.Group   = 'input';
    % Options: Scout function
    sProcess.options.scoutfunc.Comment    = {'Mean', 'Max', 'PCA', 'Std', 'All', 'Scout function:'};
    sProcess.options.scoutfunc.Type       = 'radio_line';
    sProcess.options.scoutfunc.Value      = 1;
    sProcess.options.scoutfunc.InputTypes = {'results'};
    sProcess.options.scoutfunc.Group   = 'input';
    % Options: Time-freq
    sProcess.options.edit.Comment = {'panel_timefreq_options', 'Morlet wavelet options: '};
    sProcess.options.edit.Type    = 'editpref';
    sProcess.options.edit.Value   = [];
    % Options: Normalize
    sProcess.options.labelnorm.Comment = '<BR>Spectral flattening:';
    sProcess.options.labelnorm.Type    = 'label';
    sProcess.options.normalize.Comment = {'<B>None</B>: Save non-standardized time-frequency maps', '<B>1/f compensation</B>: Multiply output values by frequency'; ...
                                          'none', 'multiply'};
    sProcess.options.normalize.Type    = 'radio_label';
    sProcess.options.normalize.Value   = 'none';
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
    
    
    
    DataToProcess = {sInputs.FileName};
    tfOPTIONS.TimeVector = in_bst(sInputs(1).FileName, 'Time');

    % === OUTPUT STUDY ===

    tfOPTIONS.iTargetStudy = [];
    
    
    
    
    % === START COMPUTATION ===
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    bin_size = 50; % This sets the bin size for the raster plot. IN MSEC % THIS SHOULD BE ON THE INITIAL GUI
    
    
    
    
    fs = abs(1./(tfOPTIONS.TimeVector(2)-tfOPTIONS.TimeVector(1)));
    
    temp = load(['E:/brainstorm_db/Playground/data/' sInputs(1).FileName],'ChannelFlag','Time');
    nElectrodes = size(temp.ChannelFlag,1); 
    raster = zeros(nElectrodes, length(sInputs), ceil(length(tfOPTIONS.TimeVector)/(bin_size*fs/1000))); % 195x3x33  = 195 electrodes x 3 trials x 33 bins.
    the_bin_centers = linspace(temp.Time(1)-abs((temp.Time(1)-temp.Time(2))/2),temp.Time(end)+abs((temp.Time(1)-temp.Time(2))/2),ceil(length(tfOPTIONS.TimeVector)/(bin_size*fs/1000)));
    clear temp
    
    for ifile = 1:length(sInputs)
        trial = load(['E:/brainstorm_db/Playground/data/' sInputs(ifile).FileName]);
        single_file_binning = zeros(nElectrodes, ceil(length(tfOPTIONS.TimeVector)/(bin_size*fs/1000)));
        for ielectrode = 1:size(trial.F,1)
            for ievent = 1:size(trial.Events,2)
                
                if strcmp(trial.Events(ievent).label, ['Spikes Electrode ' num2str(ielectrode)])
                     [~,~,bin_it_belongs_to]= histcounts(trial.Events(ievent).times,the_bin_centers);
                     a = unique(bin_it_belongs_to);
                     occurences = [a;histc(bin_it_belongs_to,a)];
                     try
                        single_file_binning(ielectrode,occurences(1,:)) = occurences(2,:);
                     catch
                         continue
                     end
                         
                end
            end
            
        end
        raster(:,ifile,:) = single_file_binning;
        
    end
    
    
    ielec = 90;
    
    figure(1);imagesc(the_bin_centers,1:length(sInputs),squeeze(raster(ielec,:,:))); 
    xlabel 'Time (sec)';ylabel 'Trial'; yticks([1:length(sInputs)]);
    title ({'Raster Plot';['Electrode: ' num2str(ielec)];['Bin Size: ' num2str(bin_size) ' msec']})
    
    
    isError = 0;
    Messages = 'Success';
    [OutputFiles] = save_virtual_freq_file(raster, the_bin_centers);
% % % % % % % %     [OutputFiles, Messages, isError] = bst_timefreq(DataToProcess, tfOPTIONS);
    if ~isempty(Messages)
        if isError
            bst_report('Error', sProcess, sInputs, Messages);
        else
            bst_report('Info', sProcess, sInputs, Messages);
            disp(['BST> process_timefreq: ' Messages]);
        end
    end
end




function [OutputFiles] = save_virtual_freq_file(raster, the_bin_centers)
    
    a = load('E:/brainstorm_db/Playground/data/Monkey/test_LFP_EYE/timefreq_trial011_morlet_170616_1605.mat'); 
    % This loads an example TF to get key information
    % (RowNames,Options,History). Get rid of that
    
    
    TF = permute(raster,[1 3 2]);
    Time = the_bin_centers;
    TFmask = true(size(raster,2),size(raster,3));
    Freqs = 1:size(TF,3);
    Std = [];
    Comment = 'Raster Plot';
    DataType = 'data';
    TimeBands = [];
    RefRowNames = [];
    RowNames = a.RowNames;
    Measure = 'power';
    Method = 'morlet';
    DataFile = 'Monkey/test_LFP_EYE/data_event1_trial011_02.mat';
    SurfaceFile = [];
    GridLoc = [];
    GridAtlas = [];
    Atlas = [];
    HeadModelFile = [];
    HeadModelType = [];
    nAvg = [];
    ColormapType = [];
    DisplayUnits = [];
    Options = a.Options;
    History = a.History;
    
    
    FileMat.TF = TF;
    FileMat.Time = Time;
    FileMat.TFmask = TFmask;
    FileMat.Freqs = Freqs;
    FileMat.Std = Std;
    FileMat.Comment = Comment;
    FileMat.DataType = DataType;
    FileMat.TimeBands = TimeBands;
    FileMat.RefRowNames = RefRowNames;
    FileMat.RowNames = RowNames;
    FileMat.Measure = Measure;
    FileMat.Method = Method;
    FileMat.DataFile = DataFile;
    FileMat.SurfaceFile = SurfaceFile;
    FileMat.GridLoc = GridLoc;
    FileMat.GridAtlas = GridAtlas;
    FileMat.Atlas = Atlas;
    FileMat.HeadModelFile = HeadModelFile;
    FileMat.HeadModelType = HeadModelType;
    FileMat.nAvg = nAvg;
    FileMat.ColormapType = ColormapType;
    FileMat.DisplayUnits = DisplayUnits;
    FileMat.Options = Options ;
    FileMat.History = History;
    
    Filename = 'E:/brainstorm_db/Playground/data/Monkey/test_LFP_EYE/timefreq_trial011_morlet_990616_1605.mat';
    OutputFiles = {Filename};
    save (Filename,'TF', 'Time', 'TFmask','Freqs','Std','Comment','DataType','TimeBands','RefRowNames',...
        'RowNames','Measure','Method','DataFile','SurfaceFile','GridLoc','GridAtlas','Atlas','HeadModelFile',...
        'HeadModelType','nAvg','ColormapType','DisplayUnits','Options','History')

%     db_add_data(iTargetStudy, FileName, FileMat);
    db_add_data(16, Filename, FileMat);

end









