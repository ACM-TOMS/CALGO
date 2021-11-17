echo on
clear
% ----------------------------------------------------------------------
% Example 9: Tensor times tensor - outer product
% ----------------------------------------------------------------------

% Create a tensor of size 3 x 4
A = tensor(reshape([1:12],[3,4]))

% Create a tensor of size 2
B = tensor([1;2],2)

% Example of tensor outer product
ttt(A,B)

echo off