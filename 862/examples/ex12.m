echo on
clear
% ----------------------------------------------------------------------
% Example 12: Matricizing a tensor using the tensor_as_matrix class.
% ----------------------------------------------------------------------

% Create a 4 x 3 x 2 tensor with the elements 1 thru 24.
T = tensor(reshape(1:24,[4 3 2]))

% Matricize T by putting its first dimension into the rows
tensor_as_matrix(T,[1])

% Matricize T by putting its second dimension into the rows
tensor_as_matrix(T,[2])

% Matricize T by putting its second dimension into the rows and
% specifying the order of the dimensions spanning the columns
tensor_as_matrix(T,[2],[3 1])

echo off
