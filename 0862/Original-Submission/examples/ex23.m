echo on
% ----------------------------------------------------------------------
% Example 23: Example of the higher-order orthogonal iteration to
% compute a Tucker approximation with a 3 x 4 x 2 core array
% ----------------------------------------------------------------------
% (Note: Must run ex20 immediately before this script.)

% Compute a Tucker approximation with a core array of size 3 x 4 x 2
% using HOOI.
T4 = hooi(T,[3 4 2])

% Convert to a dense tensor
T4f = full(T4)

% How good is the answer?
diff3 = norm(T - T4f) / norm(T)

echo off
