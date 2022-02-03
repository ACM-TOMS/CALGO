echo on
clear
% ----------------------------------------------------------------------
% Example 1: Creating a tensor object from a multidimensional array
% ----------------------------------------------------------------------

% Create a three-dimensional MDA
A = rand(3,4,2)

% Convert A to a tensor
T = tensor(A)

echo off
