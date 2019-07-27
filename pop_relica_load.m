% pop_relica_load() - load RELICA structure from disk.
%
% Usage:
%   >> pop_relica_load( EEG ); % pop up interactive window
%   >> pop_relica_load( EEG, file);
%
% Graphic interface:
%   Manually select the RELICA.mat file from disk
%
% Inputs:
%   EEG     - Input dataset
%   file    - Complete path of RELICA.mat file
%
% Author:  Fiorenzo Artoni, The Biorobotics Institute / EPFL, 2017 %
%
% References:
%
% (1) Artoni, F., Menicucci, D., Delorme, A., Makeig, S., & Micera, S. (2014).
% RELICA: a method for estimating the reliability of independent components.
% NeuroImage, 103, 391-400.

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


function [EEG, com] = pop_relica_load(EEG,file);
if nargin < 1
	help pop_relica_load;
	return;
end;
com = '';

if nargin <2
    [uifilename, uifilepath] = uigetfile('RELICA.mat','Select Relica file');
    if uifilename == 0
        disp('Cancelled by user');
        return
    end
    file = [uifilepath uifilename];
end
disp('Loading RELICA...')
a = load(file); 
EEG.etc.RELICA = a.RELICA;
disp('Done!')



