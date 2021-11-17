function A = cellDivide(M,d)
% 
% A = cellDivide(M,d)
% 
%   for M a cell of matrices the function returns the cell with each matrix
%   divided by d.
% 

m = length(M);
A = cell(1,m);

for imat = 1:m
    A{imat} = M{imat}/d;
end

% Note : speed compared with cellfun and a for loop is actually better, see
% testCellfun_for !

end