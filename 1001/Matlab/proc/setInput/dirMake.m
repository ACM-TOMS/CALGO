%% dirMake
% Set and make directories for output.
%
% Warning: This function creates a directory and a file.
%
%% Syntax
%
%   seti = dirMake(seti,usevaralpha,usevarbeta,usevardelta,usevartol)
%
%% Description
% |seti = dirMake(seti,usevaralpha,usevarbeta,usevardelta,usevartol)|
%
% * Make directory for output (seti.dirOutput).
% * Set directory name: |seti.dirname =
% 'dirOutput/<dirDatetime>_<inseti><dirSuffix><dirSuffixAdd>'|.
% * Make subfolder in output for current computation.
% * Copy setting parameters (from folder inseti) into subfolder.
%
%% Input Arguments
%
% seti : structural array
%
% _The following fields can be set by user (otherwise set default values):_
%
% * seti.gscale  : grid scaling on (1) or off (0), see <setGridScale.html>
% * seti.expData : e.g. 'fresnel' if data from Institute Fresnel are used, see <setData.html>
%
% See <setInput.html> for the following fields (they can be defined by user
% too):
%
% * seti.dirOutput 
% * seti.dirDatetime
% * seti.dirSuffix      
% * seti.dirSuffixAdd
% * seti.dirname
%
%% Output Arguments
%
% Several fields in seti, if they was not set as input argument; see above.
%
%% See Also
%
% * <setInput.html>
%
%% Code
%
function seti = dirMake(seti,usevaralpha,usevarbeta,usevardelta,usevartol)

%%
% *dirOutput*

seti = checkfield(seti,'dirOutput','output');

if exist(seti.dirOutput,'dir')~=7
    fprintf('Make dir: %s.\n',seti.dirOutput)
    mkdir(seti.dirOutput);
end

%%
% *dirDatetime*
%
% Format of date and time is yyyymmddTHHMMSS, e.g. 20160815T140912
% Note: Maybe it was set before in |varalpha|, |varbeta| or |vardelta|.
%
seti = checkfield(seti,'dirDatetime',datestr(now,'yyyymmddTHHMMSS'));

%%
% *inseti*
%
% inseti is seti.inseti

%%
% *dirSuffix*

seti = checkfield(seti,'dirSuffix','');
% (Idea: a prefix could be used too.)

%%
% *Additional Suffixes (dirSuffixAdd)*
%
% Additional suffixes are set automatically and added to dirname for
% special cases.
%
seti.dirSuffixAdd = '';

%%
% *Various input*
%
% * varalpha
% * varbeta
% * vardelta

if usevaralpha == 1
    seti.dirSuffixAdd = strcat(seti.dirSuffixAdd,'_varalpha');
end
if usevarbeta == 1
    seti.dirSuffixAdd = strcat(seti.dirSuffixAdd,'_varbeta');
end
if usevardelta == 1
    seti.dirSuffixAdd = strcat(seti.dirSuffixAdd,'_delta');
end
if usevartol == 1
    seti.dirSuffixAdd = strcat(seti.dirSuffixAdd,'_tol');
end

%%
% *Label of special cases*
%
% * grid scaling (seti.gscale)
% * use of experimental data (seti.expData), e.g. Fresnel or Simonetti
%

if isfield(seti,'gscale') && seti.gscale == 1
    seti.dirSuffixAdd = strcat(seti.dirSuffixAdd,'_gscale');
end

if isfield(seti,'expData')
    if strcmp(seti.expData,'fresnel')
        seti.dirSuffixAdd = strcat(seti.dirSuffixAdd,'_fresnel'); 
    elseif strcmp(seti.expData,'simonetti')
        seti.dirSuffixAdd = strcat(seti.dirSuffixAdd,'_simonetti'); 
    else
       error('seti.expData not found... ?!? (Neither Fresnel nor Simonetti.)');
    end
end

%%
% *Set the final dirname*

if isempty(seti.inseti)
    insetiDirName = sprintf('_noinseti');
else
    insetiDirName = sprintf('_%s',seti.inseti);
end
dirnameDefault = sprintf('%s/%s%s%s%s',...
    seti.dirOutput,seti.dirDatetime,insetiDirName,seti.dirSuffix,seti.dirSuffixAdd);
clear insetiDirName;
seti = checkfield(seti,'dirname',dirnameDefault);
clear dirnameDefault;

%%
% *Make directory for output data of current computation*

if exist(seti.dirname,'dir')~=7
    fprintf('Make dir: %s.\n',seti.dirname)
    mkdir(seti.dirname);
end

%%
% *Copy used inseti data in output folder*

if ~isempty(seti.inseti)
    src = which(seti.inseti); % source
    des = sprintf('%s/inseti_%s.m',seti.dirname,seti.inseti); % destination
    copyfile(src,des);
    clear src des;
end

end