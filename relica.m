% relica() - Estimate the reliability of independent components
%            by performing  ICA decomposition several times
%            on different data segments resampled with repetition.
%
% Usage:
%   >> EEG = relica(EEG, M,algo,mode_relica,folder_relica);
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
%   The options below require 'nsgflag' set to [1] 
%   'jobid'           - String with the client job id. This was assigned to the
%                       job when created. Use with command line option option 'run'. 
%                       Default: Prefix 'relicansg_' trailed by five digit random number. e.g 'relicansg_616402'
%   'runtime'         - Time (in hours) to allocate for running the job in NSG. 
%                       Maximun time allocation is 48 hrs. Use with command line 
%                       option option 'run'. Default: 0.5 h
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
% Beamica is part of C. Kothe's  open source BCILAB toolbox 
%
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

function [EEG]=relica(EEG,M,algo,mode_relica, folder_relica,varargin)

if isstruct(EEG)
    data = EEG.data;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% PART 1: ESTIMATION %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin <5; folder_relica = fullfile(pwd,'relicaoutput'); end
    if nargin <4; mode_relica = 'point'; end
    if nargin <3; algo = 'beamica'; end
    if nargin <2; M = 50; end
    if size(data,3)>1; mode_relica = 'point';end
    
    try
        options = varargin;
        if ~isempty( varargin )
            for i = 1:2:numel(options)
                g.(options{i}) = options{i+1};
            end
        else
            g = [];
        end
    catch
        disp('RELICA_nsg() error: calling convention {''key'', value, ... } error'); return;
    end
    try g.nsgflag;                          catch, g.nsgflag       =  0;                                               end 
    try g.icaopt;                           catch, g.icaopt        =  {};                                              end 
    if  g.nsgflag
        try g.jobid = ['relicansg_' g.jobid];   catch, g.jobid         =  ['relicansg_' num2str(floor(rand(1)*1000000))];  end
        try g.runtime;                          catch, g.runtime       =  0.5;                                             end
    else
        c = parcluster;
        try g.parpools;                         catch, g.parpools      =  c.NumWorkers;                                    end
    end
    mkdir(folder_relica); % Ramon, consider creating this folder later. In case of early failure it will be creating folders unnecessarely    
    
    %%%%%%% 1a: setup %%%%%%%%%%%%%%%%%%%%%%
    X = data(:,:);
    ltrial = size(data,2);
    sR=icassoStruct(X);
    
    %%%%%%% Check mode %%%%%%%%%%%%%%%%%%%%%
    sR.mode='both';
    %-------------------- NSG --------------------
    if ~ g.nsgflag
        
        %%%%%%% estimation %%%%%%%%%%%%%%%%%%%%%
        k=0; index=[]; tempi = zeros(1,M);
        A = cell(1,M); W = A; index = A;
        
        % Parfor stuff.
        % Check for parpool and adjust the number of workers to the value in parpools.
        hpool = gcp('nocreate');
        if isempty(hpool)
            parpool(g.parpools)
        elseif hpool.NumWorkers ~= g.parpools
            delete(gcp('nocreate'))
            parpool(g.parpools);
        end
        
        parfor i=1:M
            
            if i == 1
                X_ = X;
            else
                if strcmp(mode_relica,'trial') && size(data,3)>1
                    X_=relica_bootstrap(X,ltrial);
                else
                    X_=relica_bootstrap(X);
                end
            end
            
            in = X_(:,:);
            switch algo
                case 'beamica'
                    if i==1; nrun = 1500; else; nrun = 1200; end % ensure real ICA is good
                    [Wf,Sf] = beamica(in,{},eye(size(in,1)),pinv(sqrtm(cov(in'))),mean(in,2),nrun,0.5,0.5,false,true,true,1);
                    W_ = Wf * Sf;
                    A_ = pinv(W_);
                    
                    % Add new algorithms here
                case 'runica'
                    [Wf,Sf] = runica(in,'verbose', 'on', g.icaopt{:});
                    W_ = Wf * Sf;
                    A_ = pinv(W_);
                case 'picard'
                    if  exist('picard.m', 'file')
                        [tmp, W_] = picard(in,'verbose', true, g.icaopt{:});
                        A_ = pinv(W_);
                    else
                        error('PICARD plugin must be installed');
                    end
            end
            n=size(A_,2);
            index{i}=[repmat(i,n,1), [1:n]']';
            A{i}=A_;
            W{i}=W_;
        end
        
        % Reasigning values to sctructure sR
        sR.whiteningMatrix   = eye(size(data,1));
        sR.dewhiteningMatrix = eye(size(data,1));
        sR.index = cell2mat(index)';
        sR.A =A;
        sR.W =W;
        
        % saving file
        save([folder_relica filesep 'sR'],'sR','-v7.3')
        
    else
        % Creating structure of inputs to be saved
        % This structure is retreived later on whith the results as parameters
        % here are needed
        relicainput.eegfilename = EEG.filename;
        relicainput.M = M;
        relicainput.mode_relica = mode_relica;
        relicainput.algo = algo;
        relicainput.folder_relica = folder_relica;
        relicainput.opts = g;
        
        try
            nsg_info;  % get information on where to create the temporary file
        catch
            error('Plugin nsgportal needs to be in the MATLAB path');
        end
        
        %  Section 1: Create temporary folder and save data
        tmpJobPath = fullfile(outputfolder, 'relicatmp');
        if exist(tmpJobPath,'dir'), rmdir(tmpJobPath,'s'); end
        mkdir(tmpJobPath);
        
        % Save data in temporary folder previously created.
        % Here you may change the file name to match the one in the script you will run in NSG
        pop_saveset(EEG,'filename', EEG.filename, 'filepath', tmpJobPath);
        
        % Copy toolbox to folder. temporary until updated in NSG
        relicafolder = which('relica.m');
        relicapath = fileparts(relicafolder);
        copyfile(relicapath,tmpJobPath)
        
        % PICARD copy
        if strcmp(algo, 'picard')
            if exist('picard.m', 'file')
                picardfolder = which('picard.m');
                picardpath = fileparts(picardfolder);
                copyfile(picardpath,tmpJobPath);
            else
                error('PICARD plugin must be installed');
            end
        end
        
        % Save structure
        save(fullfile(tmpJobPath,'relicainput'),'relicainput');
        
        % Section 2
        %  Manage m-file to be executed in NSG
        % Write m-file to be run in NSG.
        % Options defined in plugin are written into the file
        
        relicansg_writejobfile
        
        % Section 3
        % Submit job to NSG
        
        pop_nsg('run',tmpJobPath,'filename', 'relicansg_job.m', 'jobid', g.jobid,'runtime', g.runtime);
        display([char(10) 'RELICA job (jobID:'  g.jobid ') has been submitted to NSG' char(10) ...
                          'Copy or keep in mind the jobID assigned to this job to retreive the results later on.' char(10)...
                          'You may follow the status of your job through pop_nsg'...
                char(10)]);
        rmdir(tmpJobPath,'s');
        return;
    end
else
    try
        nsg_info;  % get information on where to create the temporary file
    catch
        error('Plugin nsgportal needs to be in the MATLAB path');
    end
    
    pop_nsg('output',EEG);
    if exist(outputfolder,'dir')
        tmpJobPath = fullfile(outputfolder,['nsgresults_' EEG]);
        load(fullfile(tmpJobPath,'relicatmp', 'sR.mat'));
        load(fullfile(tmpJobPath,'relicatmp', 'relicainput.mat'));
        
        % Reinstating variables
        EEG  = pop_loadset(fullfile(tmpJobPath,'relicatmp',relicainput.eegfilename));
        M = relicainput.M;
        mode_relica = relicainput.mode_relica;
        algo = relicainput.algo;
        folder_relica = relicainput.folder_relica;
        g = relicainput.opts;
    else
        display('Job not finded. It may be that the job has not finished yet');
        exit;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% PART 2: PROCESSING %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M=icassoGet(sR,'M');
rdim=icassoGet(sR,'rdim');
Wt = cell2mat(sR.W');
similarity = abs(corr(Wt'));
clusterparameters={'simfcn',similarity,'s2d','sim2dis','strategy','AL','L','rdim'};
num_of_args=length(clusterparameters);

%%%%%%%%%%%%% arguments %%%%%%%%%%%%%%%%
for i=1:2:num_of_args
  switch lower(clusterparameters{i})
   case 'simfcn'
    simfcn=clusterparameters{i+1};
    % Explicit similarity matrix?
    if isnumeric(simfcn)
      if size(simfcn,1)==M && size(simfcn,2)==M
        sR.cluster.similarity=simfcn;
        sR.cluster.simfcn='<similarities given explicitly>';
      else 
        error('Explicitly given similarity matrix has wrong size!');
      end
    else
      % should be a string
      switch lower(simfcn)
       case 'abscorr'
        sR.cluster.simfcn=lower(simfcn);
       otherwise
        error('''simfcn'' must be string ''abscorr'' or an MxM similarity matrix');
      end
    end
   case 's2d'
    s2dfcn=lower(clusterparameters{i+1});
    if ~ischar(s2dfcn)
      error('''s2d'' must be a string (name of a function)');
    end
    sR.cluster.s2d=s2dfcn;
   case 'l'
    L=clusterparameters{i+1};
    if isnumeric(L)
      % The user has specified max number for clusters     
      % Check L 
      if fix(L)~=L
        error('''L'' must be an integer.');
      elseif L<2
        error('''L'' must be at least 2.');
      elseif L>M
        error('''L'' cannot be more than the number of estimates.');
      end
    else
      if ~strcmp(lower(L),'rdim')
        error('''L'' expects an integer value or ''rdim''.');
      end
      % set (reduced) data dimension
      L=icassoGet(sR,'rdim');
    end   
    if L>100
      warning(['R-index requested for more that 100 clusters: this can' ...
               ' be heavy...']);
    end
   case 'strategy'
    strategy=clusterparameters{i+1};
    if ~ischar(strategy)
      error('''strategy'' must be a string');
    end
    % we are case insensitive
    strategy=upper(strategy);
    sR.cluster.strategy=strategy;
    switch sR.cluster.strategy
     case {'AL','CL','SL'}
     % hierarchical clustering
     otherwise
      error(['Strategy ' strategy ' not implemented.']);
    end
   otherwise
    error(['Indentifier ' clusterparameters{i} ' not recognized.']);
  end
end

%%%%%%%%% Compute similarities %%%%%%%%%%
switch lower(sR.cluster.simfcn)
 case '<similarities given explicitly>'
   % already handled
 case 'abscorr'
  sR.cluster.similarity=abs(corrw(icassoGet(sR,'W'),icassoGet(sR,'dewhitemat')));
  %just to make sure  
  sR.cluster.similarity(sR.cluster.similarity>1)=1; 
  sR.cluster.similarity(sR.cluster.similarity<0)=0;
end
%%%%% Convert to dissimilarities using .s2d
D=feval(sR.cluster.s2d, sR.cluster.similarity);
%%%%% Make partition %%%%%%%%%%%%%%%
[sR.cluster.partition,sR.cluster.dendrogram.Z,sR.cluster.dendrogram.order]=...
    hcluster(D,sR.cluster.strategy);
%%%%% Compute cluster validity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% init R
sR.cluster.index.R=ones(M,1)*NaN;
% compute
sR.cluster.index.R(1:L,1)=rindex(D,sR.cluster.partition(1:L,:));  

save([folder_relica filesep 'sR'],'sR','-v7.3')
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% PART 3: PROJECT %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sR=icassoProjection(sR,'cca','s2d','sqrtsim2dis','epochs',75);
outputDimension=2;
method='cca';
projectionparameters={'s2d','sqrtsim2dis','epochs',75 ,'alpha',0.7,'epochs',75,...
	   'radius',max(icassoGet(sR,'M')/20,10),'s2d','sqrtsim2dis'};
num_of_args=length(projectionparameters);
for i=1:2:num_of_args
  switch lower(projectionparameters{i})
   case 's2d'
    sim2dis=projectionparameters{i+1};
   case 'epochs'
    epochs=projectionparameters{i+1};
   case 'alpha'
    alpha=projectionparameters{i+1};
   case 'radius'
    CCAradius=projectionparameters{i+1};
   otherwise
    error(['Indentifier ' projectionparameters{i} ' not recognized.']);
  end
end
D=feval(sim2dis,sR.cluster.similarity);
disp([char(13) 'Projection, using ' upper(method) char(13)]);
switch method 
 case 'mmds'
  P=mmds(D); 
  P=P(:,1:outputDimension);
 otherwise
  % Start from MMDS
  initialProjection=mmds(D); initialProjection=initialProjection(:,1:2);
  % 
  dummy=rand(size(D,1),outputDimension);
  % rand. init projection: set 
  % initialProjection=dummy;
  switch method
   case 'sammon' % Use SOM Toolbox Sammon 
    P=sammon(dummy,initialProjection,epochs,'steps',alpha,D);
   case 'cca'    % Use SOM Toolbox CCA 
    P=cca(dummy,initialProjection,epochs,D,alpha,CCAradius);
  end
end
sR.projection.method=method;
sR.projection.parameters=projectionparameters;
sR.projection.coordinates=P;

save([folder_relica filesep 'sR'],'sR','-v7.3')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% PART 4: POST PROCESS %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L=icassoGet(sR,'rdim');[Iq, A_centroid, W_centroid]=icassoResult(sR,L);
ncomp = size(sR.A{1},2);
A_boot_percomp{ncomp}=[]; W_boot_percomp{ncomp} =[]; indici_boot_percomp{ncomp} = [];
A_boot = sR.A; W_boot = sR.W;
nrun = length(A_boot);
indici = sR.cluster.partition(ncomp,:);
for i = 1 : length(indici)/ncomp
    b = indici(ncomp*(i-1)+1 : ncomp*i);
    for j=1:length(b) % nel cluster j ci vanno tutte le componenti 
        A_boot_percomp{b(j)}(:,end+1) = A_boot{i}(:,j);
        W_boot_percomp{b(j)}(end+1,:) = W_boot{i}(j,:);
        indici_boot_percomp{b(j)}(end+1,:) = [i j];
    end
end
indici = sR.cluster.partition(ncomp,:);
indici_boot_permatrice = reshape(indici,ncomp,length(sR.A));

real.A = sR.A{1};
real.W = sR.W{1};
real.S = sR.whiteningMatrix;
real.clustok_ord = 0;
indice = indici_boot_permatrice(:,1);
for i = 1:ncomp
    ind = find(indice == i);
    if ~isempty(ind)
        real.A_ord{i} = real.A(:,ind);
        real.Ind_ord{i} = ind;
        real.W_ord{i} = real.W(ind,:);
    else
        real.A_ord{i} = [];
        real.Ind_ord{i} = [];
        real.W_ord{i} = [];
    end        
end
RELICA.Iq = Iq;
RELICA.sR = sR;
RELICA.A_boot_percomp = A_boot_percomp;
RELICA.W_boot_percomp = W_boot_percomp;
RELICA.indici_boot_percomp = indici_boot_percomp;
RELICA.indici_boot_permatrice = indici_boot_permatrice;
RELICA.A_centroid = A_centroid;
RELICA.W_centroid = W_centroid;
RELICA.A_real = real.A;
RELICA.W_real = real.W;
RELICA.ind_real = indici_boot_permatrice(:,1); %which cluster they belong to
RELICA.folder_output = [folder_relica filesep 'RELICA'];
save([folder_relica filesep 'RELICA'],'RELICA','-v7.3');
try delete([folder_relica filesep 'sR.mat']); end
EEG.etc.RELICA = RELICA;