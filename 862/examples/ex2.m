echo on
clear
% ----------------------------------------------------------------------
% Example 2: The tensor class explicitly tracks trailing singleton
% dimensions
% ----------------------------------------------------------------------

% Create an MDA of size 4 x 3 x 1.
A = rand([4 3 1]);

% Trailing singleton dimensions are ignored, so
% the number of reported dimensions is only 2.
ndims(A)

size(A)

% Explicitly specifying the dimensions with the tensor constructor 
% creates an order-3 tensor of size 4 x 3 x 1.
T = tensor(A,[4 3 1]);

% Now the number of dimensions and size include the trailing
% singleton dimension.
ndims(T)

size(T)

% The WHOS command reports the correct sizes for the tensor
% object 
whos

echo off
