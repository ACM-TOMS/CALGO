echo on

% ----------------------------------------------------------------------
% Example 6: An alternate approach to calculating n-mode products
% ----------------------------------------------------------------------

% (Note: Must run ex5 immediately before this script.)

% Put U and V into a cell array
W{1} = U;
W{2} = V;

% Now compute A x_1 U x_3 V using the cell array
F = ttm(A,W,[1 3])

% The answers are the same
norm(C - F)

echo off
