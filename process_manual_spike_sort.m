function varargout = process_manual_spike_sort( varargin )

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
    sProcess.Comment     = 'MEA MANUAL Spike Sorter';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Exctract';
    sProcess.Index       = 1300;
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
    
    
    % Convert all the files in input
    for i = 1:length(sInputs)
        [path_to_load, filebase, ~] = bst_fileparts(file_fullpath(sInputs(i).FileName));

        wave_clus(path_to_load, filebase)
        
        
    end
end

