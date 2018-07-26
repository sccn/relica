% pop_RELICA_main() - Runs RELICA
%
% Description: The RELICA toolbox enables to estimate the reliability of
% independent components by performing the ICA algorithm several times,
% each time on different data, resampled with repetition. The algorithm is
% already optimized for GPU/CUDA, if available on the workstation.    
%
% Usage:
%   >> EEG = pop_RELICA_main(EEG); % pop up interactive window
%
% Graphic interface:
%   "Number of Bootstraps" - [edit box] enter the number of bootstrap runs,
%                   i.e. the number of times the ICA algorithm should run
%   "ICA algorithm" - [select box] select the ICA algorithm to run.
%                   Possible options are:
%                   BeamICA: Infomax-based ICA algorithm already
%                           optimized for GPU computing (if available)
%   "Select Trial by Trial or Point by Point" - [select box] select the
%                   type of bootstrap to do. Possible options are:
%                   Trial-by-trial: each time the ICA algorithm is 
%                           run trials are resampled with repetition. This                              
%                           option is only available on epoched datasets.
%                           If the dataset is no epoched, the RELICA
%                           algorithm will revert back to the
%                           point-by-point bootstrapping version
%                   Point-by-point: each time the ICA algorithm is 
%                           run data samples are resampled with repetition. 
%                           This is the default option
%   "Select Output Backup Folder" - [select box] select the output folder
%
%
% Inputs:
%   EEG     - Input dataset
%
% Outputs:
%   EEG     - Output dataset: RELICA data is in EEG.etc.RELICA
%
% Author:  Fiorenzo Artoni, The Biorobotics Institute / EPFL, 2017 %
%
% References:
%
% (1) Artoni, F., Menicucci, D., Delorme, A., Makeig, S., & Micera, S. (2014).
% RELICA: a method for estimating the reliability of independent components.
% NeuroImage, 103, 391-400.
% 

% Copyright (C) 2017 Fiorenzo Artoni, The Biorobotics Institute , EPFL, SCCN, fiorenzo.artoni@epfl.ch
%
% Clustering and relative visualization within RELICA makes use of  modified 
% routines from J. Himberg's open source FastICA - ICASSO package
% Beamica is part of C. Kothe's  open source BCILAB toolbox 
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function [EEG, com] = pop_RELICA_main(EEG);
com = '';

if nargin < 1
	help pop_RELICA_main;
	return;
end;	

userInput = inputgui('title', 'RELICA', 'geom', ...
   {{2 7 [0 0] [1 1]}   {2 7 [1 0] [1 1]} ...
    {2 7 [0 2] [1 1]}   {2 7 [1 2] [1 1]} ...
    {2 7 [0 4] [1 1]}   {2 7 [1 4] [1 1]} ...
    {2 7 [0 6] [1 1]}   {2 7 [1 6] [1 1]} ...
    {6 7 [0 9] [1 1]}},...
'uilist',...
   {{'style' 'text' 'string' 'Number of Bootstraps'} {'style' 'edit' 'string' '50','tooltipstring', wordwrap('Number of times the ICA algorithm is repeated on bootstrapped data',80)} ...
    {'style' 'text' 'string' 'ICA algorithm'}    {'style' 'popupmenu' 'string' 'beamica'}...
    {'style' 'text' 'string' 'Select Trial-By-Trial or Point-by-Point bootstrap'}    {'style' 'popupmenu' 'string' 'point-by-point|trial-by-trial'}...
    {'style' 'text' 'string' 'Select Output backup folder'}   {'style' 'edit' 'string' [pwd filesep 'RELICA'],'tooltipstring', wordwrap('Folder for saving data',80)} ...
    {'style', 'pushbutton', 'string', 'Help', 'callback', 'pophelp(''pop_RELICA_main'');'}});

if isempty(userInput)
    disp('Operation terminated by user.')
    return;
end

algo_cell = {'beamica'};
mode_relica_cell = {'point','trial'};
M = str2num(userInput{1,1}); %#ok<*ST2NM>
algo = algo_cell{userInput{1,2}};
mode_relica  = mode_relica_cell{userInput{1,3}};
folder_relica  = userInput{1,4};

[EEG] = RELICA_main(EEG,M,algo,mode_relica,folder_relica);
com = ['EEG = RELICA_main(EEG, ' num2str(M) ', ''' algo ''', ''' mode_relica ''', ''' folder_relica '''); ' ];
EEG = eegh(com, EEG);
disp('Done.')



function outtext = wordwrap(intext,nChars)
outtext = '';    
while ~isempty(intext)
    if length(intext) > nChars
        cutoff = nChars+find([intext(nChars:end) ' ']==' ',1)-1;
        outtext = [outtext intext(1:cutoff-1) '\n']; %#ok<*AGROW>
        intext = intext(cutoff+1:end);
    else 
        outtext = [outtext intext];
        intext = '';
    end
end
outtext = sprintf(outtext);


