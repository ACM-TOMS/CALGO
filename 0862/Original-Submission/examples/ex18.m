echo on
clear
% ----------------------------------------------------------------------
% Example 18: Adding two CP tensors
% ----------------------------------------------------------------------

% Create a CP tensor with one factor
A = cp_tensor(5, [2 3 4]', [1 2]', [5 4 3]');

% Add the two CP tensors, and the result has two factors
B = A + A

% Convert to a dense tensor
C = full(B)

echo off