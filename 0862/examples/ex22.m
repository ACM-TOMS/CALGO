echo on
% ----------------------------------------------------------------------
% Example 22: Example of the higher-order orthogonal iteration to
% compute a Tucker approximation with a 2 x 2 x 1 core array
% ----------------------------------------------------------------------
% (Note: Must run ex20 immediately before this script.)

% Compute a Tucker approximation with a core array of size 2 x 2 x 1
% using HOOI.
T3 = hooi(T,[2 2 1])

% Convert to a dense tensor
T3f = full(T3)

% How good is the answer?
diff3 = norm(T - T3f) / norm(T)

echo off
