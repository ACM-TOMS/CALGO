%% checkfield
% Check existence of a field in a structural array
%
%% Syntax
%
%   struct = checkfield(struct,fieldname,default)
%   struct = checkfield(struct,fieldname,default,dispDepth)
%
%% Description
% |struct = checkfield(struct,fieldname,default)|
% checks the existence of a field in the structure array |struct|, e.g. |seti| 
% (and set a default value when not).
%
% |struct = checkfield(struct,fieldname,default,dispDepth)| does the same
% as |struct = checkfield(struct,fieldname,default)| in case of |dispDepth|
% greater or equal 1, but suppresses the output message in case of
% |dispDepth = 0|.
%
%% Example
% If |seti.testdirname| does not exist we will set the parameter to something like
% |output/20160816T101411| with the following example.
%
%   seti.dirOutput = 'output';
%   seti.dirDatetime = datestr(now, 'yyyymmddTHHMMSS');
%
%   default = sprintf('%s/%s',seti.dirOutput,seti.dirDatetime);
%   seti = checkfield(seti,'testdirname',default);
%
% We will delete the fields |dirOutput| and |dirDatetime| in the following
% example
%
%   fields = {'dirOutput','dirDatetime'};
%   seti = rmfield(seti,fields);
%
%% Input Arguments
%
% * |struct|    : structure array with data, e.g. |seti| in our program
% * |fieldname| : name of the field, i.e. call is seti.(fieldname)
% * |default|   : default value (is set if seti.(fieldname) does not exist)
%
%% Output Arguments
%
% * |struct|      : structure array with data, but with struct.(fieldname)
%                 (if it was not existing before)
%
%% Code

function struct = checkfield(struct,fieldname,default,varargin)
%%

if nargin == 4
    dispDepth = varargin{1};
else
    dispDepth = 0;
end

if ~isfield(struct,fieldname)
    struct.(fieldname) = default;
    if dispDepth >= 1
        fprintf('   Set parameter %s = %s.\n',fieldname,num2str(struct.(fieldname)));
    end
end

%% From long to short version
% Long old code:
%
%   if ~isfield(seti,'dirname')
%       seti.dirname = sprintf('%s/%s_%s%s%s',seti.dirOutput,seti.dirDatetime,seti.inseti,seti.dirSuffix,seti.dirSuffixAdd);
%       message = sprintf('Parameter "dirname" was not set. Set dirname = %s.',seti.dirname);
%       disp(message)
%   end
%
% New short code:
%
%   default = sprintf('%s/%s_%s%s%s',seti.dirOutput,seti.dirDatetime,seti.inseti,seti.dirSuffix,seti.dirSuffixAdd);
%   seti = checkfield(seti,'dirname',default);
%
end