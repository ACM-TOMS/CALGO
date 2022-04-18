echo on
clear
% ----------------------------------------------------------------------
% Example 7: Comparison of ttv and ttm
% ----------------------------------------------------------------------

% Create a 3 x 4 x 2 tensor
A = tensor(rand(3,4,2));

% Create a length 3 vector
u = rand(3,1);

% Computing "tensor times vector" (ttv) yields a result of 
% size 4 x 2.
C = ttv(A,u,1)

% Computing "tensor times matrix" (ttm) yields same result,
% but of size 1 x 4 x 2.
B = ttm(A,(u'),1)

echo off