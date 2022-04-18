function M = buildProduct(A,prod)
%
% M = BUILDPRODUCT(A,PROD)
% 
% Constructs the product, of matrices in the cell A,
% corresponding to the sequence in PROD
%
% M = A{prod(end)}*A{prod(end-1)}...*A{prod(1)}
% 
% if prod(i)=0, A{prod(i)} is replaced by the identity
%


n = size(A{1},1);
l = length(prod);

M = eye(n);

for t=1:l
    if prod(t)>0
        M = A{prod(t)}*M; % The order of the multiplication is of course important and is related with the lifting: A'XA or AXA'
    end                   % in conitope for instance
end
end