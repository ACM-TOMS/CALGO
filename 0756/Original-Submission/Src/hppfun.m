function F = hppfun(y,fdat)
%HPPFUN (not intended for calling directly by the user)
%	Returns residual for solution of nonlinear equations. 
%	Used by HPPARAM.
%
%	Written by Toby Driscoll.  Last updated 5/23/95.


n = fdat(1,1);
beta = fdat(1:n-1,2);
nmlen = fdat(1:n-3,3);
rows = 1:fdat(2,1);
left = fdat(rows,4);
right = fdat(rows,5);
cmplx = fdat(rows,6);
qdat = fdat(1:fdat(3,1),7:fdat(4,1));

% Transform y (unconstr. vars) to x (prevertices)
x = [-1;cumsum([0;exp(y)])];

% Check crowding of singularities.
if any(diff(x)<eps) | any(isinf(x))
  % Since abs(y) is large, use it as the penalty function.
  F = y;
  disp('Warning: Severe crowding')
  return
end

% Compute the integrals appearing in nonlinear eqns.
xleft = x(left);
xright = x(right);
mid = mean([xleft.' ; xright.']).';
% For integrals between non-adjacent singularities, choose intermediate
% points in the upper half-plane.
mid(cmplx) = mid(cmplx) + i*(xright(cmplx)-xleft(cmplx))/2;
ints = hpquad(xleft,mid,left,x,beta,qdat) - ...
    hpquad(xright,mid,right,x,beta,qdat);

if any(ints==0)
  % Singularities were too crowded in practice.
  F = y;
  disp('Warning: Severe crowding')
else
  % Compute nonlinear equation residual values.
  F1 = abs(ints(~cmplx));		% F1(1) = abs(ints(1))
  F1 = F1(2:length(F1))/F1(1);
  F2 = ints(cmplx)/ints(1);
  F = [F1;real(F2);imag(F2)] - nmlen;
end


