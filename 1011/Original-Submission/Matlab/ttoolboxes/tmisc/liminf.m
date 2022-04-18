function v = liminf(v);
% [ w ] = liminf( v ) 
% Computes the cumulative liminf of vectors.
% w = arrayfun(@(x) min(v(x:end)),1:length(v));
%
% E.g.: liminf([10 1 9 2 8 3 7 4 6])
% 
% See also: limsup
%
% Written by: tommsch, 2018

v=arrayfun(@(x) min(v(x:end)),1:length(v));

end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 