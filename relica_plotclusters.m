% relica_plotclusters() - Plo
% Usage:
%   >> EEG = relica(EEG, M,algo,mode_relica,folder_relica);
%
% Inputs:
%   sR         - Structure sR from EEG.ETC.RELICA.sR
%
% Optional inputs:
%   nclust           - [Integer|[]] Number of clusters to partition the data. An
%                      empty array will indicate to uses as many clusters as original ICA
%                      components on the data into. Default: [] (Same number of 
%                      clusteras as original ICA components)
% labelalpha         - [0:1]. Transparency level of the background of the cluster
%                      labels. [0]: Totally transparent, [1] not transparency.
%                      Default: 0.4
% 
% labelbackcolor     - Color of the axis background. Value must be provided as
%                      a 1by3 RGB array. Default(white): [1 1 1]
% centmarkersize     - [Integer] Size of centroid marker. Default: 9
% markersize         - [Integer] Size of points. Default: 5
% axbackgroundcolor  - Background color of the figure. Value must be provided as
%                      a 1by3 RGB array. Default(light gray): [0.95 0.95 0.95]
% labelfontsize     -  Font size of the cluster labels. Default: 14

% Outputs:
%   
%
% Author:  Ramon Martinez-Cancino, UCSD, INC, SCCN 2019
%          Scott Makeig, UCSD, INC, SCCN 2019
%          Fiorenzo Artoni, 2019 
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

function hfig = relica_plotclusters(sR, varargin)

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
    disp('relica_plotclusters() error: calling convention {''key'', value, ... } error'); return;
end

try g.nclust;             catch, g.nclust             = [];                end
try g.labelalpha;         catch, g.labelalpha         = 0.4;               end
try g.labelbackcolor;     catch, g.labelbackcolor     = [1 1 1];           end
try g.centmarkersize;     catch, g.centmarkersize     = 9;                 end
try g.markersize;         catch, g.markersize         = 5;                 end
try g.axbackgroundcolor;  catch, g.axbackgroundcolor  = [0.95 0.95 0.95];  end 
try g.labelfontsize;      catch, g.labelfontsize      = 14;                end

% Number of clusters
if isempty(g.nclust)
    level=icassoGet(sR,'rdim');
else
    level = g.nclust;
end

% Check if cluster level is valid
maxCluster=size(sR.cluster.partition,1);
if level<=0 | level>maxCluster
  error('Cluster level out of range or not specified.');
end

% Get cluster membership and centroids
partition=sR.cluster.partition(level,:);
Ncluster=max(partition);
index2centrotypes=icassoIdx2Centrotype(sR,'partition',partition);

%Get projection coordinates
p=sR.projection.coordinates; 

% Visualization
%% Create colors for clusters
my_col = jet(length(unique(partition)));
my_col = my_col(randperm(length(unique(partition)), length(unique(partition))),:);
ptscolor = my_col(partition,:);

% Plotting figure
hfig = figure('Units', 'Normalized', 'Color', BACKCOLOR); hold on;
title(sprintf('Estimate space as a 2D %s projection',upper(sR.projection.method)),'Fontsize', 20);

% Axis stuff
ax = gca;
set(ax,'box','on', 'XTickLabel', '', 'YTickLabel','', 'color',g.axbackgroundcolor);
grid on;
 
for i =1:size(p,1)
    plot(p(i,1),p(i,2),'o','MarkerFaceColor',ptscolor(i,:),'MarkerEdgeColor',ptscolor(i,:), 'Markersize', g.markersize);
end

% Centering the cloud
ptsxlim = minmax(p(:,1)'); 
ptsylim = minmax(p(:,2)');
xmargin = (ptsxlim(2)-ptsxlim(1))*0.1;
ymargin = (ptsylim(2)-ptsylim(1))*0.1;

xlim([ptsxlim(1)-xmargin, ptsxlim(2)+xmargin]);
ylim([ptsylim(1)-ymargin, ptsylim(2)+ymargin]);

% Converting positions relative to figure
oldunits = get(ax, 'Units');
set(ax, 'Units', 'Normalized');
axpos = get(ax, 'Position');
set(ax, 'Units', oldunits);

% Get axes drawing area in data units
ax_xlim = xlim(ax);
ax_ylim = ylim(ax);

ax_pixels_per_xdata = axpos(3) ./ diff(ax_xlim);
ax_pixels_per_ydata = axpos(4) ./ diff(ax_ylim);

% These are figure-relative
xlabelcentroid = ((p(:,1) - ax_xlim(1)) .* ax_pixels_per_xdata + axpos(1))+1/500;
ylabelcentroid = ((p(:,2) - ax_ylim(1)) .* ax_pixels_per_ydata + axpos(2))+1/500;

% Cluster labels
txt=cellstr(num2str([1:Ncluster]'));

%% Plot centrotypes
for i=1:Ncluster
       textcoord = double([xlabelcentroid(index2centrotypes(i)), ylabelcentroid(index2centrotypes(i)), 0.02,0.02]);
       hannotation = annotation('textbox',textcoord, 'String', txt{i},'FitBoxToText', 'on', 'Margin', 0 , 'Units', 'normalized');
       set(hannotation,'HorizontalAlignment', 'center',...
                       'VerticalAlignment', 'middle',...
                       'Color', [0 0 0],...
                       'BackgroundColor', g.labelbackcolor,...
                       'FaceAlpha', g.labelalpha,...
                       'EdgeColor', 'none',...
                       'FontSize', g.labelfontsize )
                   
       plot(p(index2centrotypes(i),1),p(index2centrotypes(i),2),'-o', 'markersize', g.centmarkersize,'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor',ptscolor(index2centrotypes(i),:));
      % text(p(index2centrotypes(i),1)-xwidth/100,p(index2centrotypes(i),2),txt{i},'horiz','right','color',[0 0 0],'fontsize',15,'fontweight','bold');
end