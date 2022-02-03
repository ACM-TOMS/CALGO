function CVEC = amgkfdc(K, XBAR, X)
% AMGKFDC: CVEC = amgkfdc(K, XBAR, X)
%
% DESCRIPTION: Coefficients for finite difference approximation
% of the K'th derivative at XBAR given arbitrary locations X.
% Usually the elements X(i) are monotonically increasing and
% X(1) <= XBAR <= X(n), but neither condition is required.
% The X values need not be equally spaced but must be distinct.  
%
% INPUT: (K, XBAR, X) are required and not validated.
% K = nonnegative integer order of derivate >= 0
% XBAR = location where derivative is to be approximated
% X = abscissa locations such that NX = length(X) > K
%
% OUTPUT: CVEC = vector of coefficients such that length(CVEC) = NX.
%
% This file is part of AMGKQ.  See AMGKQ.M for details.
% Copyright (C) 2014 Robert W. Johnson, Alphawave Research.
% This is free software; see GPLv3 or greater for copying conditions.
% There is ABSOLUTELY NO WARRANTY; not even for MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE.  For details, see GPLv3 or greater.
%
% Based on FDCOEFFF.M (C) 2007 Randall J. LeVeque
% From <http://www.amath.washington.edu/~rjl/fdmbook/>
% As described in "Finite Difference Methods for Ordinary and Partial
% Differential Equations: Steady-State and Time-Dependent Problems",
% Society for Industrial and Applied Mathematics (SIAM), Philadelphia,
% Softcover / ISBN 978-0-898716-29-0, xiv+339 pages, July, 2007.
%
% Based on the program "weights" in 
%   B. Fornberg, "Calculation of weights in finite difference formulas",
%   SIAM Review 40 (1998), pp. 685-691.
%
% Note: Forberg's algorithm can be used to simultaneously compute the
% coefficients for derivatives of order 0, 1, ..., m where m <= n-1.
% This gives a coefficient matrix C(1:n,1:m) whose K'th column gives
% the coefficients for the K'th derivative.  It can be easily modified
% to return the whole array if desired.

% Consistency check
n = length(X);
if (n <= K), error('amgkfdc: n <= K'), end
m = K;
% change to m=n-1 if you want to compute coefficients for all
% possible derivatives.  Then modify to output all of C.
c1 = 1;
c4 = X(1) - XBAR;
C = zeros(n-1,m+1);
C(1,1) = 1;
for i = 1:n-1
    i1 = i+1;
    mn = min(i,m);
    c2 = 1;
    c5 = c4;
    c4 = X(i1) - XBAR;
    for j = 0:i-1
        j1 = j+1;
        c3 = X(i1) - X(j1);
        c2 = c2*c3;
        if j == i-1
            for s = mn:-1:1
                s1 = s+1;
                C(i1,s1) = c1*(s*C(i1-1,s1-1) - c5*C(i1-1,s1))/c2;
            end
        C(i1,1) = -c1*c5*C(i1-1,1)/c2;
        end
        for s = mn:-1:1
            s1 = s+1;
            C(j1,s1) = (c4*C(j1,s1) - s*C(j1,s1-1))/c3;
        end
        C(j1,1) = c4*C(j1,1)/c3;
    end
c1 = c2;
end
% last column of C gives desired coefficients
CVEC = C(:,end);  
end % function % amgkfdc

