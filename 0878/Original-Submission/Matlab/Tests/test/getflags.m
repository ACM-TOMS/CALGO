%GETFLAGS  Get flags from varargin
%
%  [ARGS, FLG1, FLG2...] = GETFLAGS(ARGS, 'FLAG1', 'FLAG2'...) sets FLG1 to true
%  if FLAG1 is a string member of the cell array ARGS, otherwise to false, and
%  similarly for FLAG2, FLAG3 etc. On exit all the FLAGS that are present have
%  been removed from ARGS.
%
%  To obtain flags with arguments (as in "plot(x,y,'linewidth',3)") use
%  [ARGS,FLG,VALUE] = GETFLAGS(ARGS, 'FLAG'). If FLAG is present in ARGS, FLG is
%  set to true and the next argument is returned in VALUE. Both FLAG and VALUE
%  are removed from ARGS. If FLAG is last in ARGS or not there, VALUE becomes
%  empty. To process several flag-value pairs use multiple calls to GETFLAGS.
%
%  EXAMPLE:
%     function somefunc(varargin)
%     [varargin, aflag, bflag] = getflags(varargin, '-a', '-b')
%     [varargin, cflag, cval] = getflags(varargin, '-c')
%
%     Assume now that somefunc is called with: somefunc('ABC','-b',2,'-c',5).
%     After the call to getflags, varargin will be {'ABC', 2}, aflag is false,
%     bflag is true, cflag is true and cval is 5. Similarly, after the call:
%                            somefunc ABC -b 2 -c 5
%     varargin will be {'ABC', '2'}, cval will be '5' and the other variables
%     will be as before.

function [args,varargout] = getflags(args,varargin)
  ascertain(nargout==nargin || nargout == 3 && nargin ==2);
  n = length(args);
  IC = false(n,1);
  for i=1:n, IC(i) = ischar(args{i}); end;
  Keep = true(n,1);
  Jc = find(IC);
  for k=1:length(varargin)
    Iq = strmatch(varargin{k},args(IC),'exact');  
    Keep(Jc(Iq)) = false;
    varargout{k} = ~isempty(Iq);
    if nargout > nargin
      if ~isempty(Iq) && length(args) > Jc(Iq),
        varargout{2} = args{Jc(Iq)+1};
        Keep(Jc(Iq)+1) = false;
      else
        varargout{2} = [];
      end
    end
  end
  args = args(Keep);
end
