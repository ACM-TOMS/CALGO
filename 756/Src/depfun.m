function F = depfun(y,fdat)
%DEPFUN (not intended for calling directly by the user)
%	Returns residual for solution of nonlinear equations.
%	Used by DEPARAM.
%
%	Written by Toby Driscoll.  Last updated 5/23/95.

n = fdat(1,1);
beta = fdat(1:n,2);
nmlen = fdat(1:n-3,3);
qdat = fdat(1:fdat(2,1),4:fdat(3,1));

% Transform y (unconstr. vars) to z (prevertices)
cs = cumsum(cumprod([1;exp(-y)]));
theta = 2*pi*cs(1:n-1)/cs(n);
z = ones(n,1);
z(1:n-1) = exp(i*theta);

% Check crowding.
if any(diff(theta)<eps) | any(isnan(theta))
  % Since abs(y) is large, use it as the penalty function.
  F = y;
  disp('Warning: Severe crowding')
  return
end

% Compute the integrals appearing in nonlinear eqns.
dtheta = theta(2:n-1) - theta(1:n-2);
dtheta(dtheta > pi) = dtheta(dtheta > pi) - 2*pi;
mid = exp(i*(theta(1:n-2) + dtheta/2));

ints = dequad(z(1:n-2),mid,1:n-2,z,beta,qdat) - ...
    dequad(z(2:n-1),mid,2:n-1,z,beta,qdat);

if any(ints==0)
  % Singularities were too crowded in practice.
  F = y;
  disp('Warning: Severe crowding')
else
  % Compute equation residual values.
  F = abs(ints(2:n-2))/abs(ints(1)) - nmlen;  

  % Compute residue.
  res = -sum(beta./z)/ints(1);

  F = [F;real(res);imag(res)];
end


