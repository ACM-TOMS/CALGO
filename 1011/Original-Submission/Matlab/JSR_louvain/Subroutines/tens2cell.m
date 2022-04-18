function M = tens2cell(A)
%
% TENS2CELL transforms 3-D tensor to cell array of matrices
% 
%  M = TENS2CELL(A)
%    returns M a cell array of matrices such that M{i} = A(:,:,i)
%         

m = size(A,3);
M = cell(1,m);

for i=1:m
    M{i} = A(:,:,i);
end

end