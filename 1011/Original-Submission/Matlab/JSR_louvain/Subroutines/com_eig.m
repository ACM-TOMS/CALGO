function [v, a, b] = com_eig(A, B, x)

% COM_EIG Computes a common eigenvector of two matrices.
%
%    [v, a, b] = com_eig(A, B) tries to find a common vector of the two matrices
%      A, B, using a Newton method. If successfull, the vector is returned in v
%      and a, b are the corresponding eigenvalues, else v is [].
%
%    [v, a, b] = com_eig(A, B, x) does the same using x as starting point for
%      the Newton method.
%
% REFERENCES
%    A.E.Ghazi, S.E.Hajji, L.Giraud, S.Gratton,
%      "Newton's method for the common eigenvector problem",
%      Journal of Computational and Applied Mathematics, 219:398-407, 2008

v = [];
n = size(A, 1);
if (nargin < 3),
    x = rand(n, 1);
end

tol = 1e-12;
x = x/norm(x);
I = eye(n);
maxiter = 80;
retry = 10;

while retry > 0,
    y = x;
    a = y'*A*x;
    b = y'*B*x;
    e = A*x - a*x;
    f = B*x - b*x;
    iter = 0;

    while ((norm(e) > tol) || (norm(f) > tol)) && (iter < maxiter),
        F = [e ; f ; y'*x - 1];
        J = [ (I - x*y')*A - a*I ;
              (I - x*y')*B - b*I ;
                     y'         ];

        x = x - J\F;
        y = x/norm(x);
        a = y'*A*x;
        b = y'*B*x;
        e = A*x - a*x;
        f = B*x - b*x;
        iter = iter + 1;
    end
    if (iter ~= maxiter),
        v = x;
        return;
    end
    retry = retry - 1;
    x = rand(n, 1);
end