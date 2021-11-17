echo on
clear
% ----------------------------------------------------------------------
% Example 5: Calculating n-mode products
% ----------------------------------------------------------------------

% Create a 4 x 3 x 2 tensor
A = tensor(rand(4,3,2));

% Create a 2 x 4 matrix
U = rand(2,4);

% Create a 3 x 2 matrix
V = rand(3,2);

% Computing A x_1 U x_3 V
B = ttm(A,U,1);
C = ttm(B,V,3)

% Computing A x_3 V x_1 U
D = ttm(A,V,3);
E = ttm(D,U,1)

% The answers are the same
norm(C - E)

echo off