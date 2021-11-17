echo on
clear
% ----------------------------------------------------------------------
% Example 3: The tensor object can explicitly store zero- and
% one-dimensional objects.
% ----------------------------------------------------------------------

% A scalar can be stored as a tensor of order zero.
T0 = tensor(5,[]);

order(T0)

ndims(T0)

size(T0)

% A vector can be stored as a tensor of order one.
T1 = tensor(rand(4,1),[4]);
order(T1)

ndims(T1)

size(T1)

% The WHOS command does not report the correct sizes!
whos

echo off
