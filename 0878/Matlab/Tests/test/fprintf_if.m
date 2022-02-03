% FPRINTF_IF  Print if condition
%
%   FPRINTF_IF(COND, FMT, P1, P2,...) prints P1, P2,... with format FMT if the
%   condition COND is true, otherwise nothing is printed.

function fprintf_if(cond, varargin)
  if cond
    fprintf(varargin{:})
  end
end