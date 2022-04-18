echo on
clear
% ----------------------------------------------------------------------
% Example 10: Computing the contracted product of a tensor times a
% tensor 
% ----------------------------------------------------------------------

% Create a 3 x 4 x 2 tensor
A = tensor(rand(3,4,2))

% Create a 4 x 3 x 2 tensor
B = tensor(rand(4,3,2))

% Tensor multiplication - matching dimensions 1 & 3 of A with
% dimensions 2 & 3 of B. 
ttt(A,B,[1 3],[2 3])

echo off
