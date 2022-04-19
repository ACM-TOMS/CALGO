function [ N ] = ndimsm(X)
% [ N ] = ndimsm(X)
% Consistent behaviour of ndims, with regards to multi-dimensional applications.
% ndimsm = length(sizem(X)).
%
% E.g.: ndimsm([1; 2; 3])
%
% See also: ndims
%
% Written by: tommsch, 2017

if(nargin==0 || isempty(X)) ; 
    N=0; 
else
    N = length(sizem(X));
end

end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 

