echo on
clear
% ----------------------------------------------------------------------
% Example 20: Higher-order power method for approximation
% ----------------------------------------------------------------------

% Create a tensor
T = tensor(rand(3,4,2))

% Compute an approximation with one factor
T1 = hopm(T)

% Convert T1 to a densor tensor
T1f = full(T1)

% How good is the approximation?
diff1 = norm(T - T1f) / norm(T)

echo off