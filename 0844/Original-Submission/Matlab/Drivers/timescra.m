% This program compares scra with the Matlab sparse svd program.
% The matrix is generated from a set of singular values using
% the matlab funciton sprandn.  For the meaning of the parameters
% see the prologue of testspqr.
%
% Coded by G. W. (Pete) Stewart
% Jun 24 2004


% Set parameters.

m = 10000;
n = 10000;
p = min(m,n);
tol = 1e-3;
svalmin = 6;
gappos = 10;
gapfac = 1e-6;

% Generate the matrix.  For large values of m and n, this is
% a time consuming process.

S = logspace(0,-svalmin,p);
%S(s:n) = gapfac*S(gappos:n);
A = sprandn(m, n, .001, S);
disp('end generation')
pause

% Time the programs.

for ncr=10:5:40

disp(ncr)

tic
[nc, cx, nr, rx, T, rsd] = scra(A, tol, ncr);
toc
disp([nc, nr])

tic
[U,s,V] = svds(A,ncr);
toc
end
