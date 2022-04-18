function mu=constructmu(k,dim)
% mu = constructmu(k, dim)
% Makes all differences up to order k, each column is one difference, or equivalently
% mmakes all multiindices of dimension dim which sum up to k
%
% Input:
%   k       order
%   dim     dimension
%
% Output:
%   mu      multiindices
%
% Eg.: constructmu(2,2)
%
% See also: differencescheme
%
% Written by: tommsch, 2016

mu=mixvector(0:k,dim);  %all possible multiindices which sum is less than k
[~,osum]=sort(sum(mu,1));
mu=mu(:,osum);
mu=mu(:,sum(mu,1)==k & sum(mu,1)>=0);
mu=fliplr(mu);

end


function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 