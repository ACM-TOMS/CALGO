% This is a test script for spqr.  Given m and n, it constructs
% an mxn matrix A with singular values 1, ..., 10^-svalmin.  It
% the singular values numbered gp, ..., p (p = min([m,n])) are multiplied
% by gapfac to create a gap.  It then runs spqr with values in the
% workspace for its arguments.
%
% Input to testspqr:
%
% m, n:     The size of the matrix
% svalmin:  Used a described above to compute the singular values of
%           A
% gappos:   The position of the gap.  What is acutally used is
%           gp = min([gappos, m, n])
% gapfac:   The factor by which singular values gp, ..., p is
%           multiplied.  Set to 1 for no gap.
%
% Coded by G. W. (Pete) Stewart
% Jun 23 2004

% Set up the matrix.

p = min([m,n]);
if m >= n
   [U, temp] = qr(randn(m,n),0);
   [V, temp] = qr(randn(n));
else
   [U, temp] = qr(randn(m));
   [V, temp] = qr(randn(n,m),0);
end
if p>1
   S = logspace(0,-svalmin,p);
else
   S = [1];
end
gp = min([gappos,p]);
S(gp:p) = gapfac*S(gp:p);
S = diag(S);
A = U*S*V';

% Print out the values used by testspqr

disp(sprintf('m, n = %3d, %3d', m , n))
disp(sprintf('svalmin, gappos, gapfac = %3d, %3d, %9.1e', ...
              svalmin, gappos, gapfac))

% Print out the arguments to spqr

disp(sprintf('tol, maxcols, fullR, pivot, cn = %9.1e, %3d, %3d, %3d, %3d', ...
              tol, maxcols, fullR, pivot, cn))

% Call spqr

[ncols, R, colx, norms] = spqr(A, tol, maxcols, fullR, pivot, cn);

% Print the results

disp(sprintf('colx'))
disp(sprintf('%3d', colx))
disp(sprintf('ncols = %3d', ncols))
if size(norms,1) ~= 0
   disp(sprintf('norms(1:ncols)'))
   disp(sprintf('%9.1e', norms(1:ncols)))
   if ncols < n
      disp(sprintf('norms(ncols+1:n)'))
      disp(sprintf('%9.1e', norms(ncols+1:n)))
   end
end

% Compute and print diagnostic statistics.

Y = A(:,colx(1:ncols));
Q = Y/R(1:ncols,1:ncols);
disp(sprintf('Orthogonality measure for Q = %9.1e', norm(eye(ncols) - Q'*Q)))
disp(sprintf('Reproduction of A (R = Q''A) = %9.1e', norm(A - Q*(Q'*A))))
if fullR
disp(sprintf('Reproduction of A (R from spqr) = %9.1e', norm(A(:,colx) - Q*R)))
end
disp(' ')
