% EEG = RELICA_main(EEG,M,algo,mode_relica,folder_relica) - RUNS RELICA
%
% Description: The RELICA toolbox enables to estimate the reliability of
% independent components by performing the ICA algorithm several times,
% each time on different data, resampled with repetition and letting you see. The algorithm is
% already optimized for GPU/CUDA, if available on the workstation.    
%
% Usage:
%   >> EEG = RELICA_main(EEG, M,algo,mode_relica,folder_relica); % pop up interactive window
%
% Inputs:
%   EEG         - Input dataset
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
% Outputs:
%   EEG     - Output dataset: RELICA data is in EEG.etc.RELICA, the same is
%             saved in folder_relica folder
%
% Author:  Fiorenzo Artoni, The Biorobotics Institute / EPFL, 2017 %
%
% Reference:
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

function [EEG]=RELICA_main(EEG,M,algo,mode_relica,folder_relica)

data = EEG.data;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% PART 1: ESTIMATION %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin <5;folder_relica = [pwd filesep 'temp']; end; mkdir(folder_relica);
if nargin <4;mode_relica = 'point'; end;
if nargin <3;algo = 'beamica'; end;
if nargin <2;M = 50; end;
if size(data,3)>1; mode_relica = 'point';end


%%%%%%% 1a: setup %%%%%%%%%%%%%%%%%%%%%%
X = data(:,:);
l_trial = size(data,2);
sR=icassoStruct(X); 

%%%%%%% Check mode %%%%%%%%%%%%%%%%%%%%%
sR.mode='both';

%%%%%%% estimation %%%%%%%%%%%%%%%%%%%%%
k=0; index=[]; tempi = zeros(1,M);
for i=1:M
  tic
  if i < 3 
    fprintf('\n\n%s\n\n',['Randomization using ' algo ' Round ' num2str(i) '/' num2str(M) '... Time estimate not availabe yet.']);
  else
    t_taken = seconds(mean(tempi(2:end))*(M-i+1)); [h,m,s] = hms(t_taken); 
    fprintf('\n\n%s\n\n',['Randomization using ' algo ' Round ' num2str(i) '/' num2str(M) '... Remaining ' num2str(h) ' Hours, ' num2str(m) ' Minutes.']);   
  end
  if i == 1
      X_ = X;
  else
      if strcmp(mode_relica,'trial') && size(data,3)>1
          X_=bootstrap_trial(X,ltrial);         
      else
          X_=bootstrap(X);
      end
  end
  sR.whiteningMatrix = eye(size(data,1));
  sR.dewhiteningMatrix = eye(size(data,1));
  switch algo
      case 'beamica'
         if i==1; nrun = 1500; else; nrun = 1200; end % ensure real ICA is good
         in = X_(:,:);
         [Wf,Sf] = beamica(in,{},eye(size(in,1)),pinv(sqrtm(cov(in'))),mean(in,2),nrun,0.5,0.5,false,true,true,1);
         W_ = Wf * Sf;
         A_ = pinv(W_);  
      % add new algorithms here   
  end
  n=size(A_,2);
  k=k+1;
  sR.index(end+1:end+n,:)=[repmat(k,n,1), [1:n]'];
  sR.A{k}=A_; sR.W{k}=W_; 
  tempi(i)=toc;
end
save([folder_relica filesep 'sR'],'sR','-v7.3')
         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% PART 2: PROCESSING %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M=icassoGet(sR,'M');
rdim=icassoGet(sR,'rdim');
clusterparameters={'simfcn','abscorr','s2d','sim2dis','strategy','AL','L','rdim'};
num_of_args=length(clusterparameters);

%%%%%%%%%%%%% arguments %%%%%%%%%%%%%%%%
for i=1:2:num_of_args;
  switch lower(clusterparameters{i})
   case 'simfcn'
    simfcn=clusterparameters{i+1};
    % Explicit similarity matrix?
    if isnumeric(simfcn),
      if size(simfcn,1)==M & size(simfcn,2)==M,
        sR.cluster.similarity=simfcn;
        sR.cluster.simfcn='<similarities given explicitly>';
      else 
        error('Explicitly given similarity matrix has wrong size!');
      end
    else
      % should be a string
      switch lower(simfcn)
       case 'abscorr'
        ; % ok
        sR.cluster.simfcn=lower(simfcn);
       otherwise
        error('''simfcn'' must be string ''abscorr'' or an MxM similarity matrix');
      end
    end
   case 's2d'
    s2dfcn=lower(clusterparameters{i+1});
    if ~ischar(s2dfcn),
      error('''s2d'' must be a string (name of a function)');
    end
    sR.cluster.s2d=s2dfcn;
   case 'l'
    L=clusterparameters{i+1};
    if isnumeric(L),
      % The user has specified max number for clusters     
      % Check L 
      if fix(L)~=L,
        error('''L'' must be an integer.');
      elseif L<2,
        error('''L'' must be at least 2.');
      elseif L>M,
        error('''L'' cannot be more than the number of estimates.');
      end
    else
      if ~strcmp(lower(L),'rdim'),
        error('''L'' expects an integer value or ''rdim''.');
      end
      % set (reduced) data dimension
      L=icassoGet(sR,'rdim');
    end   
    if L>100,
      warning(['R-index requested for more that 100 clusters: this can' ...
               ' be heavy...']);
    end
   case 'strategy'
    strategy=clusterparameters{i+1};
    if ~ischar(strategy),
      error('''strategy'' must be a string');
    end
    % we are case insensitive
    strategy=upper(strategy);
    sR.cluster.strategy=strategy;
    switch sR.cluster.strategy
     case {'AL','CL','SL'}
      ; % hierarchical clustering
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
  ; % already handled
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
for i=1:2:num_of_args;
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
try; delete([folder_relica filesep 'sR.mat']); end;
EEG.etc.RELICA = RELICA;
% if isempty(EEG.icaweights)
%     EEG.icaweights = RELICA.W_real;
%     EEG.icawinv = RELICA.A_real;
%     EEG.icasphere = eye(size(EEG.data,1));
% end

function X=bootstrap(X)
N=size(X,2);
index=round(rand(N,1)*N+.5);
X=X(:,index);

function X=bootstrap_trial(X,ltrial)
N=size(X,2)/ltrial;
T = reshape(X,size(X,1),ltrial,N);
index=round(rand(N,1)*N+.5);
T=T(:,:,index);
X = reshape(T,size(T,1),prod(size(T(1,:))));















