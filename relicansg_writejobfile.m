% Author:  Ramon Martinez-Cancino, UCSD, INC , SCCN  2019
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
% Defininng variables
fprintf(fid, 'mode_relica   = ''%s'';\n', mode_relica);
fprintf(fid, 'M             = %u;\n', M);
fprintf(fid, 'algo          = ''%s'';\n', algo);
fprintf(fid, 'folder_relica = eval(''fullfile(pwd, ''''relicaoutput'''')'');\n');
fprintf(fid, 'icaopt     = eval(''%s'');\n', opttext);

% Calling relica
fprintf(fid,'[EEG]=relica(EEG,M,algo,mode_relica, folder_relica,''local'',''icaopt'', icaopt);\n');
fprintf(fid,'pop_saveset(EEG,''filename'', EEG.filename, ''filepath'', pwd)');
fclose(fid);