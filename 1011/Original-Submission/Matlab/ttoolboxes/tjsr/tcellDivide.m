function A = tcellDivide(M,d)
% 
% A = tcellDivide(M,d)
% 
%   for M a cell of matrices the function returns the cell with each matrix
%   divided by d.
% 
% Original version: October 2018
% Modified by: tommsch, January 2019

m = length(M);
A = cell(1,m);

for imat = 1:m
    A{imat} = M{imat}/d;
end

% Note : speed compared with cellfun and a for loop is actually better, see
% testCellfun_for !

end