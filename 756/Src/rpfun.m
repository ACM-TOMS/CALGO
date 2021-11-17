function F = rpfun(y,fdat)
%RPFUN  (not intended for calling directly by the user)
%	Returns residual for solution of nonlinear equations. 
%	Used by RPARAM.
%
%	Written by Toby Driscoll.  Last updated 5/23/95.

n = fdat(1,1);
beta = fdat(1:n+2,2);
nmlen = fdat(1:fdat(2,1)-1,3);
rows = 1:fdat(2,1);
left = fdat(rows,4);
right = fdat(rows,5);
cmplx = fdat(rows,6);
qdat = fdat(1:fdat(3,1),7:fdat(4,1));
corners = fdat(5:8,1);

% Transform y (unconstr. vars) to z (actual params)
z = rptrnsfm(y,corners);
nb = corners(3)-1;

% Check crowding of singularities.
if any(abs(diff(z))<eps) | any(isinf(z))
  F = y;
  disp('Warning: Severe crowding')
  return
end

% Compute the integrals appearing in nonlinear eqns.
zleft = z(left);
zright = z(right);
mid = mean([zleft.' ; zright.']).';
ints = stquad(zleft,mid,left,z,beta,qdat) - ...
    stquad(zright,mid,right,z,beta,qdat);

if any(ints==0)|any(isnan(ints))
  % Singularities were too crowded.  Try returning a
  % big residual to get fsolve to try something else.
  %%F = 10*(5+randn(1))*ones(n-3,1);
  F = y;
  disp('Warning: Severe crowding')
else
  % Compute nonlinear equation residual values.
  cmplx2 = cmplx(2:length(cmplx));
  F1 = abs(ints(~cmplx)); 		% F1(1) = abs(ints(1))
  F1 = log( (F1(2:length(F1))/F1(1)) ./ nmlen(~cmplx2) );
  F2 = log( (ints(cmplx)/ints(1)) ./ nmlen(cmplx2) );
  F = [F1;real(F2);imag(F2)];
end


