function [ nn ] = ordering2double( nn, base)
% d = ordering2double( oo, [base] )
% Takes a vector of numbers an returns an integer having this number representation in backwards ordering
%
% Input:
%   nn      vector, vector of numbers
%   base    int, default=10, base in which the output number is represented. def
%
% Output:
%   d       the number
%
% Note: 
%   The order of the numbers gets mirrored, e.g.: ordering2double([1 2 3 4]) gives 4321.
%
% E.g.: n = ordering2double([1 5 2 3 2 2 ])



if(nargin==1); base=10; end;
nn=sum(nn(:).*repmat(base,[numel(nn),1]).^((0:numel(nn)-1)'));

end