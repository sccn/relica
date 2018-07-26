%
% Description: Plots the results of the RELICA functions    
%
% Usage:
%   >> EEG = pop_RELICA_plots(EEG); % pop up interactive window
%
% Graphic interface:
%   "Plot Graph" - [check box] plot the clusters in a 2D CCA space. The
%           larger the cluster the less stable the corresponding IC. The
%           picture does not make sense without a good number of
%           repetitions (e.g., 50)
%   "Plot real maps" - [check box] plot the scalp maps of the
%           non-bootstrapped ICA performed within RELICA. Each scalp map
%           bears the number of the cluster it belongs to               
%   "Plot cluster bootstrap maps"  - [check box] plot the bootstrapped
%           maps of a cluster.
%   "Which cluster?"  - [edit box]  number of the cluster to plot
%   "How many plots?" - [edit box] limit the number of scalp maps to plot
%
% Inputs:
%   EEG     - Input dataset
%
% Outputs:
%   EEG     - Output dataset: RELICA data is in EEG.etc.RELICA
%
% Author:  Fiorenzo Artoni, The Biorobotics Institute / EPFL, 2017 %

% Reference:
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

function [EEG, com] = pop_RELICA_plots(EEG);
if nargin < 1
	help pop_RELICA_plots;
	return;
end;

userInput = inputgui('title', 'pop_RELICA_plots()', 'geometry',  ...
   {1 1 1 [3 1] [3 1] 1},...
'uilist',...
   {{ 'style' 'checkbox' 'string' 'plot graph' 'value' 0} ...
    { 'style' 'checkbox' 'string' 'plot real maps' 'value' 0} ...
    { 'style' 'checkbox' 'string' 'plot cluster bootstrap maps' 'value' 0}...    
    {'style' 'text' 'string' 'which cluster?'}   {'style' 'edit' 'string' [1],'tooltipstring', wordwrap('Bootstrap cluster to plot',80)} ...
    {'style' 'text' 'string' 'how many plots?'}   {'style' 'edit' 'string' [20],'tooltipstring', wordwrap('Bootstrap cluster to plot',80)} ...
    {'style', 'pushbutton', 'string', 'Help', 'callback', 'pophelp(''pop_RELICA_plots'');'}});

if isempty(userInput)
    error('Operation terminated by user.')
end

plot_graph = userInput{1,1}; 
plot_real = userInput{1,2}; 
plot_clust = userInput{1,3};
cls = str2num(userInput{1,4}); 
cls_nplots = str2num(userInput{1,5});

if plot_graph
    RELICA_plots(EEG,'cluster')
    com = ['EEG = RELICA_plots(EEG, ''cluster'');' ];
    EEG = eegh(com, EEG);
end
if plot_real
    RELICA_plots(EEG,'real_maps')
    com = ['EEG = RELICA_plots(EEG, ''real_maps'');' ];
    EEG = eegh(com, EEG);
end
if plot_clust
    RELICA_plots(EEG,'cls_maps',cls,cls_nplots)
    com = ['EEG = RELICA_plots(EEG, ''cls_maps'', ' num2str(cls) ', ' num2str(cls_nplots) ');' ];
    EEG = eegh(com, EEG);
end

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


