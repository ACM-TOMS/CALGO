echo on
clear
% ----------------------------------------------------------------------
% Example 19: Creating a Tucker tensor
% ----------------------------------------------------------------------

% Create a core array
lambda = tensor(rand(4,3,1),[4 3 1]);

% Create the factors
for n = 1 : 3
    U{n} = rand(5,size(lambda,n));
end

% Assemble into a Tucker tensor
A = tucker_tensor(lambda,U)

% The size of A
size(A)

% The size of the core array
size(A.lambda)

echo off