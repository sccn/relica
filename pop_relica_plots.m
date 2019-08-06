%
% Description: Plots the results of the RELICA functions    
%
% Usage:
%   >> EEG = pop_relica_plots(EEG); % pop up interactive window
%
% Graphic interface:
% Plot type [popup menu]      - Type of plots (see options below)
%   "IC cluster grap"         - plot the clusters in a 2D CCA space. The
%                               larger the cluster the less stable the corresponding IC. The
%                               picture does not make sense without a good number of
%                               repetitions (e.g., 50)
%   "Real IC maps"            - plot the scalp maps of the
%                               non-bootstrapped ICA performed within RELICA. Each scalp map
%                               bears the number of the cluster it belongs to               
%   "Cluster bootstrap maps"  - plot the bootstrapped maps of a cluster.
%
%   "Cluster to plot"    - [edit box]  number of the cluster to plot
%   "Max number of maps" - [edit box] limit the number of scalp maps to plot
%   "Plot"               - [Button] Plot selected option
%
% Inputs:
%   EEG     - Input dataset
%
% Outputs:
%   EEG     - Output dataset: RELICA data is in EEG.etc.RELICA
%
% Author:  Fiorenzo Artoni, The Biorobotics Institute / EPFL, 2017 %
%          Ramon Martinez-Cancino, SCCN, INC, UCSD

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

function [EEG, com] = pop_relica_plots(EEG)
com = '';

if nargin < 1
	help pop_relica_plots;
	return;
end

% Callbacks
cbclose = 'close(gcbf)';

cbplot = ['plottype   =  get(findobj(''tag'', ''relicaplotlist''),''value'');'...
          'cls        = str2double(get(findobj(''tag'', ''relicaplot_cls''),''string''));'...
          'cls_nplots = str2double(get(findobj(''tag'', ''relicaplot_nmaps''),''string''));'...
          'switch plottype,'...
              'case 1, relica_plots(EEG,''cluster''); com = ''EEG = relica_plots(EEG, ''''cluster'''');'';'...
              'case 2, relica_plots(EEG,''real_maps''); com = ''EEG = relica_plots(EEG, ''''real_maps'''');'';'...
              'case 3, relica_plots(EEG,''cls_maps'',cls,cls_nplots); com = [''EEG = relica_plots(EEG, ''''cls_maps'''' , '' num2str(cls) '', '' num2str(cls_nplots) '');'' ];'...
          'end;'...
          'EEG = eegh(com, EEG);'...
          'tmp = get(gcbf,''userdata'');'...
          'set(gcbf,''userdata'', [tmp com]);'];

cbpopup = 'set(findobj(''tag'', ''relicaplot_cls''), ''enable'', ''on''); set(findobj(''tag'', ''relicaplot_nmaps''), ''enable'', ''on'');if get(findobj(''tag'', ''relicaplotlist''), ''value'')~=3, set(findobj(''tag'', ''relicaplot_cls''), ''enable'', ''off''); set(findobj(''tag'', ''relicaplot_nmaps''), ''enable'', ''off''); end';

% GUI elements
plotlist = {'Bootstrap IC clusters', 'Cluster exemplar maps' 'Cluster scalp maps'};
uilist = {{ 'style' 'text' 'string' 'Plot type'} {'style'  'popupmenu'  'string' plotlist 'tag' 'relicaplotlist' 'callback' cbpopup}...
    {'style' 'text' 'string' 'Cluster to plot'} {'style' 'edit' 'string' [1],'tooltipstring', wordwrap('Bootstrap cluster to plot',80), 'tag', 'relicaplot_cls', 'enable', 'off'}...
    {'style' 'text' 'string' 'Max number of maps'} {'style' 'edit' 'string' [20],'tooltipstring', wordwrap('Maiximun number of IC maps from cluster specified to plot',80),  'tag', 'relicaplot_nmaps',  'enable', 'off'}...
    {'style'  'pushbutton' 'string' 'Help', 'tag' 'help' 'callback', 'pophelp(''pop_relica_plots'')' }...
    {'style'  'pushbutton' 'string' 'Plot' 'tag' 'plotbutton' 'callback' cbplot} };
wt = 4; ht=4;
geom = {{wt ht [0 0] [1 1] } {wt ht [1 0] [2 1] } ...
    {wt ht [0 1] [1 1] } {wt ht [1 1] [1 1] } {wt ht [2 1] [1 1] } {wt ht [3 1] [1 1] } ...
    {wt ht [0 3] [1 1] }                      {wt ht [2 3] [2 1] }};

% Create GUI
[tmp, com] = inputgui( 'geom', geom, 'uilist', uilist, 'title' , 'Plot RELICA results -- pop_relica_plots()', 'userdata',' ');

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