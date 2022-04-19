echo on
% ----------------------------------------------------------------------
% Example 21: Example of the higher-order orthogonal iteration to
% compute a Tucker approximation with a 1 x 1 x 1 core array
% ----------------------------------------------------------------------
% (Note: Must run ex20 immediately before this script.)

% Compute a Tucker approximation with a core array of size 1 x 1 x 1
% using HOOI.
T2 = hooi(T,[1 1 1])

% Convert to a dense tensor
T2f = full(T2)

% How good is the answer?
diff2 = norm(T - T2f) / norm(T)

echo off
