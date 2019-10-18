% eegplugin_RELICA() - EEGLAB plugin for estimating the Reliability of
%                       Independent Components                  
% Usage:
%   >> eegplugin_RELICA(fig, trystrs, catchstrs);
%
% Inputs:
%   fig        - [integer]  EEGLAB figure
%   trystrs    - [struct] "try" strings for menu callbacks.
%   catchstrs  - [struct] "catch" strings for menu callbacks.
%
% Author:  Dr. Fiorenzo Artoni EPFL, 2019 %
%
% References:
% (1) Artoni, F., Menicucci, D., Delorme, A., Makeig, S., & Micera, S. (2014).
% RELICA: a method for estimating the reliability of independent components.
% NeuroImage, 103, 391-400.          
% 
% (2) Artoni, F., Delorme A., Makeig S. (2018) 
% Applying dimension reduction to EEG data by Principal Component Analysis
% reduces the quality of its subsequent Independent Component
% decomposition, Neuroimage 175 176-187
%
% This project was in part supported by the European Union's Horizon 2020
% research and innovation programme under Marie Sklodowska-Curie Action
% agreement no. 750947 (BIREHAB)
% 
% Acknowledgments go to Ramon Martinez-Cancino (SCCN/INC/UCSD 2019) for making the
% algorithm available and parallelized on the NSG server and including other ICA algorithms
% and Arnaud Delorme and Scott Makeig (SCCN/INC/UCSD 2019) for the constant inputs
% and ideas to perfect the functionality.
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

function vers = eegplugin_RELICA( fig, try_strings, catch_strings); 
vers = 'RELICA1.0';
% create menu
toolsmenu = findobj(fig, 'tag', 'tools');
submenu = uimenu( toolsmenu, 'label', 'RELICA');
% add new submenu
uimenu( submenu, 'label', 'Run RELICA', 'callback', 'EEG = pop_relica(EEG);');
uimenu( submenu, 'label', 'Load RELICA from disk', 'callback', 'EEG = pop_relica_load(EEG);');
uimenu( submenu, 'label', 'Plot results', 'callback', 'EEG = pop_relica_plots(EEG);');
uimenu( submenu, 'label', 'Apply RELICA to main dataset', 'callback', 'EEG.icawinv = EEG.etc.RELICA.A_real; EEG.icaweights = EEG.etc.RELICA.W_real; EEG.icasphere = eye(size(EEG.data,1)); EEG= eeg_checkset(EEG); disp(''Applying weights to dataset... Done!'')');
