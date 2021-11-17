% SCRA produces a reduced-rank approximation to a matrix.  Specifically,
% 
%       [nc, cx, nr, rx, T, err] = cra(A, tol, maxnrc)
% 
%    produces an approximation of the form
% 
%       A(:,cx)*T*A(rx,:)
% 
%    where T is a nc by nr matrix.  The parameter err is a bound
%    on the accuracy of the approximation and is contrlled by tol.
%    The parameter maxnrc is an upper bound on nc
%    and nr.
% 
%    SCRA uses SPQR to compute SPQR factorizations of A and A'.
%
% Author: G. W. (Pete) Stewart, Jun 24 2004

function [nc, cx, nr, rx, T, err] = cra(A, tol, maxrc)

[nc, R, cx, cn] = spqr(A, tol, maxrc, 0);
[nr, S, rx, rn] = spqr(A', tol, maxrc, 0);


cx = cx(1:nc);
rx = rx(1:nr);

T = [];
for i=1:nr
   T = [T, full(A(:,cx)'*(A*A(rx(i),:)'))];
end

T = R\((R'\(T/S))/S');


err = sqrt(cn(nc)^2 + rn(nr)^2);
