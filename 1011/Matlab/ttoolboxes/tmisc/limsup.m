function v = limsup(v);
% [ w ] = limsup( v ) 
% Computes the cumulative limsup of vectors.
% w = arrayfun(@(x) max(v(x:end)),1:length(v));
%
% E.g.: limsup([10 1 9 2 8 3 7 4 6])
% 
% See also: liminf
%
% Written by: tommsch, 2018

v=arrayfun(@(x) max(v(x:end)),1:length(v));

end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 