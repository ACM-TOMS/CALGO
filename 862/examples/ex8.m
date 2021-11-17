echo on
clear
% ----------------------------------------------------------------------
% Example 8: Tensor times a sequence of vectors
% ----------------------------------------------------------------------

% Create a 3 x 4 x 2 tensor
A = tensor(rand(3,4,2));

% Create three vectors in a cell array
U = {rand(3,1), rand(4,1), rand(2,1)};

% A tensor times a sequence of vectors in every
% dimension produces a scalar.
ttv(A,U)

% A tensor times a sequence of vectors in every
% dimension *except one* produces a vector.  Here 
% we neglect the 2nd dimension.
ttv(A,U,-2)

echo off
