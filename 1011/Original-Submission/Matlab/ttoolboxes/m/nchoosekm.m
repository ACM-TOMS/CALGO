function [ z ] = nchoosekm(x,y)
% [ y ] = nchoosekm( x, y)
% Multivariate binomial coefficient
%
% E.g.: nchoosekm([10 5],[5 2])
%           
% See also: factorialm
%
% Written by tommsch, 2018

for i=1:length(x)
    x(i) = prod(nchoosek(x(i),y(i)));
end
z=prod(x);


end