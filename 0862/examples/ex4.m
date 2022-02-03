echo on
clear
% ----------------------------------------------------------------------
% Example 4: Accessors and assignment for a tensor object work the
% same as they would for a multidimensional array.
% ----------------------------------------------------------------------

% Create a random 2 x 2 x 2 tensor
A = tensor(rand(2,2,2))

% Access the (2,1,1) element
A(2,1,1)

% Reassign a 2 x 2 submatrix to be
% the 2 x 2 identity matrix
A(:,1,:) = eye(2)

% Access the A(:,1,:) sub-tensor
A(:,1,:)

echo off