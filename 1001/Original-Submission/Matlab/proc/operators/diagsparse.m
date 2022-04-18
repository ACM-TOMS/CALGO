%% diagsparse
% Creates a sparse matrix by writing the entries of a vector on the 
% diagonal.
%
%% Syntax
%
%   diagh = diagsparse(h)
%
%% Description
% |diagh = diagsparse(h)| creates a sparse matrix, |diagh|, 
% by writing the entries of the vector |h| on the diagonal.
%
%% Example
%
%   h = [1 4 5];
%   diagh = diagsparse(h);
%
% _Result:_
%  diagh =
%
%   (1,1)        1
%   (2,2)        4
%   (3,3)        5
%
%% Input Argument
%
% * h : a vector
%
%% Output Argument
%
% * diagh : matrix (stored as sparse matrix)
%
%% Code
function diagh = diagsparse(h)

% method 1: slow
% tic; diag1 = diag(h); toc

% method 2: fast
%tic
i = 1:length(h);
diagh = sparse(i,i,h);
%toc

end