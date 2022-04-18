function F = stpfun(y,fdat)
%STPFUN (not intended for calling directly by the user)
%	Returns residual for solution of nonlinear equations. 
%	Used by STPARAM.
%
%	Written by Toby Driscoll.  Last updated 5/23/95.

n = fdat(1,1);
nb = fdat(2,1);
beta = fdat(1:n+2,2);
nmlen = fdat(1:fdat(3,1)-1,3);
rows = 1:fdat(3,1);
left = fdat(rows,4);
right = fdat(rows,5);
cmplx = fdat(rows,6);
qdat = fdat(1:fdat(4,1),7:fdat(5,1));

% In this function, n refers to the number of FINITE prevertices.

% Transform y (unconstr. vars) to z (actual params)
z = zeros(n,1);
z(2:nb) = cumsum(exp(y(1:nb-1)));
z(nb+1:n) = i+cumsum([y(nb);-exp(y(nb+1:n-1))]);

% Check crowding of singularities.
if any(abs(diff(z))<eps) | any(isinf(z))
  % Try to make fsolve take a different step.
  F = y;
  disp('Warning: Severe crowding')
  return
end

% Compute the integrals appearing in nonlinear eqns.
zleft = z(left);
zright = z(right);
mid = mean([zleft.' ; zright.']).';
c2 = cmplx;
c2(2) = 0;
mid(c2) = mid(c2) - sign(left(c2)'-nb)*i/2;
ints = stquad(zleft,mid,left,z,beta,qdat) - ...
    stquad(zright,mid,right,z,beta,qdat);

absval = abs(ints(~cmplx)); 		% absval(1) = abs(ints(1))
if ~absval(1)
  rat1 = 0;
  rat2 = 0;
else
  rat1 = absval(2:length(absval))/absval(1);
  rat2 = ints(cmplx)/ints(1);
end

if any([rat1;rat2]==0) | any(isnan([rat1;rat2])) | any(isinf([rat1;rat2]))
  % Singularities were too crowded.  Try returning a
  % big residual to get fsolve to try something else.
  F = y;
  disp('Warning: Severe crowding')
else
  % Compute nonlinear equation residual values.
  cmplx2 = cmplx(2:length(cmplx));
  F1 = log( rat1 ./ nmlen(~cmplx2) );
  F2 = log( rat2 ./ nmlen(cmplx2) );
  F = [F1;real(F2);imag(F2)];
end


