% pop_relica() - Call interactive windows of RELICA Toolbox
%                The RELICA toolbox enables to estimate the reliability of
%                independent components by performing the ICA algorithm several times,
%                each time on different data, resampled with repetition and letting 
%                you see. The algorithm is already optimized for GPU/CUDA, 
%                if available on the workstation.
%
%
% Usage:
%   >>  EEG = pop_relica(EEG,M,algo,mode_relica, folder_relica);                 % Local computation call
%   >>  EEG = pop_relica(EEG,M,algo,mode_relica, folder_relica, 'runtime', 0.5); % NSG computation call
%   >>  EEG = pop_relica(EEG);                                                   % pop up interactive window
%   >>  EEG = pop_relica('relicansg_myjobID');                                   % Retreive/continue RELICA 
%                                                                                  computation performed at NSG
%
% Inputs:
%   EEG         - Input dataset in the case of running RELICA on a local computational
%                 resource or when submitting a RELICA computaion to NSG.  
%                 If retreiving a RELICA computation from NSG, this inpust
%                 should be the job ID assigned to the job.
%   M           - Number of boostraps, i.e. the number of times the ICA algorithm should run               
%   algo        - ICA algorithm to run. Currently possible options are:
%               'beamica': (default)Infomax-based ICA algorithm already
%                        optimized for GPU computing (if available)
%   mode_relica - select thd type of bootstrap to do. Possible options are:
%               'trial': each time the ICA algorithm is 
%                        run trials are resampled with repetition. This                              
%                        option is only available on epoched datasets.
%                        If the dataset is no epoched, the RELICA
%                        algorithm will revert back to the
%                        point-by-point bootstrapping version
%                'point' Point-by-point: each time the ICA algorithm is 
%                        run data samples are resampled with repetition. 
%                        This is the default option
%   folder_relica - select the output folder.
%
% Optional inputs:
% These options must be provided as a  pair: 'optname', optvalue
%   'icaopt'          - Cell array with options for the ICA method used. The format 
%                       of the inputs must be provided as pairs in a cell array as 
%                       the  following example. { 'optname1', optvalue1,'optname2', optvalue2  }
%   'parpools'        - Number of workers to use in the parallelization. 
%                       The default is the maximum number of MATLAB workers in your
%                       system (usually the number of cores). This option
%                       can not be used if the option below ('nsgflag') is
%                       set to [1].
%   'nsgflag'         - [0|1] Flag to enable [1] or disable [0] computation
%                       on NSG. Default: 0
% The options below require 'nsgflag' set to [1] 
%   'jobid'           - String with the client job id. This was assigned to the
%                       job when created. Use with command line option option 'run'. 
%                       Default: Prefix 'relicansg_' trailed by five digit random number. e.g 'relicansg_616402'
%   'runtime'         - Time (in hours) to allocate for running the job in NSG. 
%                       Maximun time allocation is 48 hrs. Use with command line 
%                       option option 'run'. Default: 0.5
% GUI Inputs:
%    'ICA method'      - Correspond to input 'algo'
%    'Relica mode'     - Correspond to input 'mode_relica'
%                        'trial-by-trial' -> 'trial'
%                        'point-by-point' -> 'point'
%    'Bootstraps'      - Correspond to input 'M'
%    'Output folder'   - Correspond to input 'folder_relica'
%    'NSG options'     - RELICA options non specific to NSG ('jobid', 'runtime')
%    'Compute on NSG'  - Correspond to input 'nsgflag'
%    'NSG options'     - Currently allows only options 'jobid' and 'runtime'. 
%                        The number of options may increase in the feature.
% 
% 'Outputs:
%       EEG     - Output dataset: RELICA data is in EEG.etc.RELICA, the same is
%                 saved in folder_relica folder
%
% See also: relica.m
%
% Authors: Fiorenzo Artoni, The Biorobotics Institute / EPFL, 2019
%          Ramon Martinez-Cancino  SCCN/INC/UCSD 2019
          
% Copyright (C) 2019 Ramon Martinez-Cancino and Fiorenzo Artoni
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

function [EEGout, com] = pop_relica(EEG, M,algo,mode_relica, folder_relica,varargin)

if nargin < 1
   help pop_relica;
	return;
end	
com = ''; EEGout = struct;

% Checking if NSGPORTAL Toolbox is set
nsginstalled_flag = 1;
try
    nsg_info;  % From here, get information on where to create the temporary file
catch
    nsginstalled_flag = 0;
end

if isstruct(EEG)
    
    icaalgorithm  = {'beamica' 'runica' 'picard'}; % add new agorithm here
    
    % Check if PICARD exist and adding it to menu
    if exist('picard.m', 'file')
        icaalgorithm{3}  = 'picard'; 
    end
    
    moderelicaval = {'point' 'trial'};
    jobiddef      = ['relicansg_' num2str(floor(rand(1)*1000000))];
    
    if license('test', 'Distrib_Computing_Toolbox')
        clust = parcluster;
        defpools = clust.NumWorkers;
    end
    
    if ~exist('M', 'var'),             M             = 50;                            end
    if ~exist('algo', 'var'),          algo          = 'beamica';                     end
    if ~exist('mode_relica', 'var'),   mode_relica   = 'point';                       end
    if ~exist('folder_relica', 'var'), folder_relica = fullfile(pwd,'relicaoutput');  end
    
    g = finputcheck(varargin, ...
                            {'icaopt'           'cell'     ''                    {} ; ...  % NSG specific
                             'jobid'            'string'   ''              jobiddef ; ...  % NSG specific
                             'runtime'          'real'     [0.1 48]        0.5 ; ...       % NSG specific
                             'nsgflag'          'real'     [0 1]           0});            % NSG specific
    if ischar(g), error(g); end
    
    if nargin < 5
        
        % Closing open GUI and creating a new one
        openfig = findobj('tag', 'pop_relicagui');
        if ~isempty(openfig)
            disp('pop_relica warning: there can be only one pop_relica window, closing old one...')
            close(openfig);
        end
        
        guititle = 'Estimate reliability of independent components -- pop_relica()';
        cfolder = '';
        nsgcheck = ['if ' num2str(nsginstalled_flag) ',if get(findobj(''tag'',''chckbx_nsgt''),''value''),' ...
                        'set(findobj(''tag'',''nsgopt''),''enable'', ''on'');'...
                     'else, set(findobj(''tag'',''nsgopt''),''enable'', ''off'');end;'...
                   'else, set(findobj(''tag'',''nsgopt''),''enable'', ''off''); set(findobj(''tag'',''chckbx_nsgt''),''value'',0); end'];
               
        cbfolder       = 'pathname = uigetdir();if ~isequal(pathname, 0),[trash,filename] = fileparts(pathname); set(findobj(gcbf, ''tag'', ''relicafolder''), ''string'', pathname);end;';
        moderelicagui = {'point-by-point' 'trial-by-trial'};
        
        uilist = {{'style' 'text' 'string' 'ICA method'}     {'style' 'popupmenu'  'string' icaalgorithm  'tag' 'icamethod' 'value' 1} ...
            {'style' 'text' 'string' 'RELICA mode'}    {'style' 'popupmenu'  'string' moderelicagui 'tag' 'relicamode' 'value' 1}...
            {'style' 'text' 'string' 'Bootstraps'}     {'style' 'edit'       'string' num2str(M)    'tag' 'nboots'}...
            {'style' 'text' 'string' 'Output folder'}  {'style' 'edit'       'string' folder_relica 'tag' 'relicafolder'} {'style' 'pushbutton' 'string' 'Browse...'   'tag' 'filebrowse' 'callback' cbfolder}...
            {'style' 'text' 'string' 'RELICA options'} {'style' 'edit'       'string' ' '           'tag' 'relicaopt'}...
            {'style' 'text' 'string' 'Compute on NSG'} {'style' 'checkbox' 'tag' 'chckbx_nsgt' 'callback' nsgcheck 'value' g.nsgflag}...
            {'style' 'text' 'string' 'NSG options'}    {'style' 'edit'       'string' ' '           'tag' 'nsgopt' 'enable' fastif( g.nsgflag, 'on','off')}};
        %%
        ht = 6; wt = 4.5;   c1 = 0; c2 = 1; c3 = 3; c4 = 4;
        geom = { {wt ht [c1 0]  [1 1]}  {wt ht [c2 0] [3.5 1]} ...
            {wt ht [c1 1]  [1 1]}  {wt ht [c2 1] [1.3 1]} {wt ht [2.3 1] [1 1]} {wt ht [3.1 1] [1 1]} ...
            {wt ht [c1 2]  [1 1]}  {wt ht [c2 2] [2.6 1]} {wt ht [3.5 2] [1 1]}...
            {wt ht [c1 3]  [1 1]}  {wt ht [c2 3] [3.5 1]}...
            {wt ht [c1 4.5]  [1 1]}  {wt ht [c2 4.5] [1 1]} ...
            {wt ht [c1 5.5]  [1 1]}  {wt ht [c2 5.5] [3.5 1]} ...
            };
        
        result = inputgui('title', guititle, 'geom', geom, 'uilist',uilist, 'eval', 'set(gcf,''tag'', ''pop_relicagui'')', 'helpcom','pophelp(''pop_relica'');');
        if length(result) == 0 return; end
        errorflag = 0;
        algo = icaalgorithm{result{1}};
        mode_relica = moderelicaval{result{2}};
        if ~isempty(result{3}), M = str2num(result{3}); else, errorflag=1; end
        if ~isempty(result{4}) && exist(result{4},'dir'), folder_relica = result{4}; else, errorflag=1; end
         
        % NSG options
        args = {'nsgflag',result{6}};
        if result{6}
            tmpoptparams = eval(['{' result{7} '}']);
            if ~isempty(tmpoptparams)
                count = 1;
                for i =1:length(tmpoptparams)/2
                    args = { args{:},  tmpoptparams{count}, tmpoptparams{count+1}};
                    count = count + 2;
                end
            end
        end
        
        % RELICA options
        tmpoptparams = eval(['{' result{5} '}']);
        if ~isempty(tmpoptparams)
            count = 1;
            for i =1:length(tmpoptparams)/2
                if ~ismember( tmpoptparams{count},{'jobid' 'runtime' 'nnodes'})
                args = { args{:},  tmpoptparams{count}, tmpoptparams{count+1}};
                count = count + 2;
                else
                    errordlg2(['Invalid NSG option provided.' char(10) 'Only RELICA options are accepted in ''RELICA options''']);
                    error(['Invalid NSG option provided.' char(10) 'Only RELICA options are accepted in ''RELICA otpions''']);
                end
            end
        end
        
        if size(EEG.data,3)>1; mode_relica = 'point';end
        
    else % no interactive inputs
        args = varargin;
    end
    
    EEGout = relica(EEG,M,algo,mode_relica,folder_relica,args{:});
else
    if ~ nsginstalled_flag
        error('Plugin nsgportal needs to be in the MATLAB path');
    end
    EEGout = relica(EEG);
end

% Command history output
if isstruct(EEG)
    tmparg = ['''' args{1} ''',' num2str(args{2})];
    count = 3;
    for iarg = 1: length(args)/2-1
        if ~iscell(args{count+1})
            tmparg = [tmparg ',''' args{count} ''',' num2str(args{count+1})];
        else           
            opttext = '{';
            for i = 1:length(args{count+1})
                if isstr(args{count+1}{i})
                    opttext = [opttext  '''' args{count+1}{i} ''''];
                elseif isnumeric(args{count+1}{i})
                    opttext = [opttext ' ' num2str(args{count+1}{i})];
                elseif islogical(args{count+1}{i})
                    opttext = [opttext ' ' num2str(args{count+1}{i})];
                end
            end
            opttext = [opttext '}'];
            tmparg = [tmparg ',''' args{count} ''',' opttext];
        end
        count = count + 2;
    end
    com = ['EEG = pop_relica(EEG, ' num2str(M) ', ''' algo ''', ''' mode_relica ''', ''' folder_relica ''',' tmparg '); ' ];
else
    com = 'EEG = pop_relica(EEG)';
end
EEGout = eegh(com, EEGout);
disp('Done.')
end