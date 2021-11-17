function F = dpfun(y,fdat)
%DPFUN  (not intended for calling directly by the user)
%       Returns residual for solution of nonlinear equations.
%       Used by DPARAM.
%
% Written by Toby Driscoll.  Last updated 5/26/95.

n = fdat(1,1);
beta = fdat(1:n,2);
nmlen = fdat(1:n-3,3);
rows = 1:fdat(2,1);
left = fdat(rows,4);
right = fdat(rows,5);
cmplx = fdat(rows,6);
qdat = fdat(1:fdat(3,1),7:fdat(4,1));

% Convert y values to z (prevertices)
cs = cumsum(cumprod([1;exp(-y)]));
theta = pi*cs(1:n-3)/cs(length(cs));
z = ones(n,1);
z(1:n-3) = exp(i*theta);
z(n-2:n-1) = [-1;-i];

% Check crowding.
if any(diff(theta)<eps) | any(isnan(theta))
  % Since abs(y) is large, use it as the penalty function.
  F = y;
  disp('Warning: Severe crowding')
  return
end

% Compute the integrals appearing in nonlinear eqns.
zleft = z(left);
zright = z(right);
mid = mean([zleft.' ; zright.']).';
% For integrals between nonadjacent singularities, choose 0 as intermediate
% integration point.
mid(cmplx) = zeros(size(mid(cmplx)));
ints = ( dquad(zleft,mid,left,z,beta,qdat) - ...
    dquad(zright,mid,right,z,beta,qdat) );

if any(ints==0)
  % Singularities were too crowded in practice.
  F = y;
  disp('Warning: Severe crowding')
else
  % Compute nonlinear equation residual values.
  F1 = abs(ints(~cmplx)); 		% F1(1) = abs(ints(1))
  F1 = F1(2:length(F1))/F1(1);
  F2 = ints(cmplx)/ints(1);
  F = [F1;real(F2);imag(F2)] - nmlen;
end

