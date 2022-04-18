function [ oo, pp, nrep ] = vector2ordering(varargin)
% [ oo, pp, nrep ] = vector2ordering(oo, [options])
% Wrapper function for findperiod.
%
% See also: findperiod
%
% Written by: tommsch, 2018

[ oo, pp, nrep ]=findperiod(varargin{:});

end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 