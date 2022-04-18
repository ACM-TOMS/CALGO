function M = tbuildproduct_fast(A,oo,v)
% M = tbuildproduct_fast(A, oo)
% Constructs the product of matrices in the cell A, corresponding to the sequence in prod.
% M = A{prod(end)}*A{prod(end-1)}...*A{prod(1)}
%
% Input:
%   A       cell array of matrices
%   oo      ordering in which to multiply them
%   v       optional. If given the function computes M*v
%           this may be faster and more accurate
%
% Output:
%   M       the product
%
% Info:
%   if prod(i)=0, A{prod(i)} is replaced by the identity
%
% E.g.: tbuildproduct_fast({2 4 5},[1 0 3])
%
% See also: tbuildproduct
% 
% Written by: Jungers, Taken from: JSR-Toolbox


n = size(A{1},1);
l = length(oo);

if(nargin==3)
    M=v;
else
    M = eye(n);
end

for t=1:l
    if oo(t)>0
        M = A{oo(t)}*M; % The order of the multiplication is of course important and is related with the lifting: A'XA or AXA'
    end                   % in conitope for instance
end
end