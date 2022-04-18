echo on
clear
% ----------------------------------------------------------------------
% Example 17: Creating a CP tensor
% ----------------------------------------------------------------------

% Create a CP tensor of size 3 x 2 x 3
A = cp_tensor(5, [2 3 4]', [1 2]', [5 4 3]')

% The size of A
size(A)

% The number of factors of A
length(A.lambda)

% Convert to a (dense) tensor
B = full(A)

echo off