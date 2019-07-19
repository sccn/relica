
% Testing vars
% tmpJobPath = pwd;
% mode_relica =  'point';
% M = 200;
% parpools = 4;
% algo = 'beamica';

%-----------

opttext = '{';
for i = 1:length(g.icaopt)
    if isstr(g.icaopt{i})
        opttext = [opttext  ' ''''' g.icaopt{i} ''''''];
    elseif isnumeric(g.icaopt{i})
        opttext = [opttext ' ' num2str(g.icaopt{i})];
    elseif islogical(g.icaopt{i})
        opttext = [opttext ' ' num2str(g.icaopt{i})];
    end
end
opttext = [opttext '}'];

fid = fopen( fullfile(tmpJobPath,'relicansg_job.m'), 'w');
fprintf(fid, 'eeglab;\n');
fprintf(fid, 'EEG  = pop_loadset(''%s'');\n', EEG.filename);
fprintf(fid, 'data = EEG.data;\n');
fprintf(fid, 'X    = data(:,:);\n');
fprintf(fid, '\n');
fprintf(fid, 'mode_relica = ''%s'';\n', mode_relica);
fprintf(fid, 'M           = %u;\n', M);
fprintf(fid, 'algo        = ''%s'';\n', algo);
fprintf(fid, 'icaopt     = eval(''%s'');\n', opttext);
fprintf(fid, '\n');
fprintf(fid, 'ltrial      = size(data,2);\n');
fprintf(fid, 'sR = icassoStruct(X);\n');
fprintf(fid, 'sR.mode =''both'';\n');
fprintf(fid, '\n');
fprintf(fid, 'k=0; index=[]; tempi = zeros(1,M);\n');
fprintf(fid, 'A = cell(1,M); W = A; index = A;\n');
fprintf(fid, '\n');
fprintf(fid, 'clusterinfo = parcluster;\n');
fprintf(fid, 'parpool(clusterinfo.NumWorkers);\n');
fprintf(fid, 'parfor i=1:M\n');
fprintf(fid, 'if i == 1\n');
fprintf(fid, 'X_ = X;\n');
fprintf(fid, 'else\n');
fprintf(fid, 'if strcmp(mode_relica,''trial'') && size(data,3)>1 \n');
fprintf(fid, 'X_=relica_bootstrap(X,ltrial);\n');
fprintf(fid, 'else\n');
fprintf(fid, 'X_=relica_bootstrap(X);\n');
fprintf(fid, 'end\n');
fprintf(fid, 'end\n');
fprintf(fid, 'in = X_(:,:);\n');
fprintf(fid, 'switch algo \n');
% BEAMICA
fprintf(fid, 'case ''beamica''\n');
fprintf(fid, 'if i==1; nrun = 1500; else; nrun = 1200; end\n');
fprintf(fid, '[Wf,Sf] = beamica(in,{},eye(size(in,1)),pinv(sqrtm(cov(in''))),mean(in,2),nrun,0.5,0.5,false,true,true,1);\n');
fprintf(fid, 'W_ = Wf * Sf;\n');
fprintf(fid, 'A_ = pinv(W_);\n');
fprintf(fid, '\n');
% RUNICA
fprintf(fid, 'case ''runica''\n');
fprintf(fid, '[Wf,Sf] = runica(in,''verbose'', ''on'', icaopt{:});');
fprintf(fid, 'W_ = Wf * Sf;\n');
fprintf(fid, 'A_ = pinv(W_);\n');
fprintf(fid, '\n');
% PICARD
fprintf(fid, 'case ''picard''\n');
fprintf(fid, 'if  exist(''picard.m'', ''file'')\n');
fprintf(fid, '[tmp, W_] = picard(in,''verbose'', true, icaopt{:});\n');
fprintf(fid, 'A_ = pinv(W_);\n');
fprintf(fid, 'end\n');

fprintf(fid, 'end\n');
fprintf(fid, 'n=size(A_,2);\n');
fprintf(fid, 'index{i}=[repmat(i,n,1), [1:n]'']'';\n');
fprintf(fid, 'A{i}=A_; W{i}=W_;\n');
fprintf(fid, 'end\n');
fprintf(fid, '\n');
fprintf(fid, 'sR.whiteningMatrix   = eye(size(data,1));\n');
fprintf(fid, 'sR.dewhiteningMatrix = eye(size(data,1));\n');
fprintf(fid, 'sR.index = cell2mat(index)'';\n');
fprintf(fid, 'sR.A =A; sR.W =W;\n');
fprintf(fid, '\n');
fprintf(fid, 'save(''sR'',''sR'',''-v7.3'')\n');
fclose(fid);