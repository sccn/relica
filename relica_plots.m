% EEG = relica_plots(EEG,graphtype,cls,cls_nplots)
%
% Description: Plots the results of the RELICA functions    
%
% Usage:
%   >> EEG = relica_plots(EEG,graphtype,cls,cls_nplots);
%
% Inputs:
%   EEG         - Input dataset
%   graphtype   - Type of graph to plot:
%                  'cluster': plot the clusters in a 2D CCA space. The
%                       larger the cluster the less stable the corresponding IC. The
%                       picture does not make sense without a good number of
%                       repetitions (e.g., 50)
%                  'real_maps': plot the scalp maps of the
%                       non-bootstrapped ICA performed within RELICA. Each scalp map
%                       bears the number of the cluster it belongs to
%                  'cls_maps': plot the bootstrapped
%                       maps of a cluster.
%   cls         - Only if cls_maps is selected: plot the boostrapped maps
%                       of a cluster
%   cls_nplots  - limit the number of scalp maps to plot
% Optional inputs:
%   sortsi      - [0,1] Sort maps by the estimated stability index. [0] Do
%                  not sort, [1] sort maps. Default: 0
%
% Outputs:
%   EEG     - Output dataset: RELICA data is in EEG.etc.RELICA, the same is
%             saved in folder_relica folder
%
% Author:  Dr. Fiorenzo Artoni, 2019 %
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
% algorithm available and parallelized on the NSG server.
% Acknowledgments go also to Arnaud Delorme and Scott Makeig (SCCN/INC/UCSD 2019) 
% for the precious inputs and ideas to perfect the project.
% Clustering and relative visualization within RELICA makes use of  modified 
% routines from J. Himberg's open source FastICA - ICASSO package
% Beamica is part of C. Kothe's  open source BCILAB toolbox %
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


function EEG = relica_plots(EEG,graphtype,cls,cls_nplots, varargin)

if nargin<3;cls = [];cls_nplots = [];end
RELICA = EEG.etc.RELICA;
icadefs;

try
    options = varargin;
    if ~isempty( varargin )
        for i = 1:2:numel(options)
            g.(options{i}) = options{i+1};
        end
    else, g= [];
    end
catch
    disp('relica_plots() error: calling convention {''key'', value, ... } error'); return;
end
try g.sortsi;   catch, g.sortsi   = 0;   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%     CLUSTERS       %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(graphtype,'cluster')
    figure;
    hcluster = relica_plotclusters(EEG.etc.RELICA.sR); % icassoGraph(RELICA.sR,'line','off','hull','off');
    set(hcluster,'name','RELICA: clusters','Color',BACKEEGLABCOLOR);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%     MAPS           %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(graphtype,'real_maps')
    if isempty(EEG.chanlocs)
        error('RELICA requires channel locations to plot the scalp maps!')
    end
    n_figs = ceil(size(RELICA.A_real,2)/20);
    if g.sortsi
        [~,iqindx] = sort(RELICA.Iq, 'descend');
    end
    for i =1:n_figs
        figure;
        set(gcf,'name',['RELICA: real maps fig' num2str(i) '/' num2str(n_figs)]);
        subp = 0;
        
        if g.sortsi
            comporder = iqindx;
        else
            comporder = 1:size(RELICA.A_real,2);
        end
        for j =20*i-19 :min(20*i,size(RELICA.A_real,2))  
            subp = subp + 1;
            n_cls = RELICA.ind_real(comporder(j));
            if sum(RELICA.ind_real==n_cls)>1;quality = 'm';else;quality = '';end
            subplot(4,5,subp);topoplot(RELICA.A_real(:,comporder(j)),EEG.chanlocs,'electrodes','off');
            title([' Cls ' num2str(n_cls) ' (' num2str(round(RELICA.Iq(comporder(j))*100)) '%)'  ]);
           %title([num2str(j) ' Cls-' num2str(n_cls) ' (' num2str(round(RELICA.Iq(n_cls)*100)) '%)'  ])
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%     MAPS  CLS      %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(graphtype,'cls_maps')
    if isempty(EEG.chanlocs)
        error('RELICA requires channel locations to plot the scalp maps!')
    end
    if ~isempty(cls)
        A = RELICA.A_boot_percomp{cls};
        if isnan(cls_nplots);nplot = size(A,2);else;nplot = min(size(A,2),cls_nplots);end        
        n_figs = ceil(nplot/20);
        for i =1:n_figs
            figure;
            set(gcf,'name',['RELICA: maps cluster' num2str(cls) ' fig ' num2str(i) '/'  num2str(n_figs)]);
            subp = 0;
            for j =20*i-19 :min(20*i,nplot)  
                subp = subp + 1;                
                subplot(4,5,subp);topoplot(A(:,j),EEG.chanlocs,'electrodes','off');
                title(['Cls ' num2str(cls) ' IC ' num2str(j)])
            end
        end
    end
end     
 