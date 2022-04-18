function [yp,yprime] = r2strip(zp,z,L)
%R2STRIP Map from rectangle to strip.
%	R2STRIP(ZP,Z,L) maps from a rectangle to the strip 0 <= Im z
%	<= 1, with the function log(sn(z|m))/pi, where sn is a Jacobi
%	elliptic function and m = exp(-2*pi*L).  The prevertices of the
%	map (in the rectangle domain) are given by Z; only the corners
%	of the rectangle defined by Z are used.
%
%       The functionality is NOT parallel to HP2DISK and DISK2HP.
%
%	Written by Toby Driscoll.  Last updated 5/23/95.

% Uses different forms based on conformal modulus of the rectangle to avoid
% underflow when modulus is large (as measured by L, separation between
% corner images on the strip).  Also returns the derivative of the map at
% the given points.

K = max(real(z));
Kp = max(imag(z));
yp = zp;
yprime = zp;
if L < 5.9
  m = exp(-2*pi*L);
  [sn1,cn1,dn1] = ellipj(real(zp),m);
  [sn2,cn2,dn2] = ellipj(imag(zp),1-m);
  sn = (sn1.*dn2 + i*sn2.*cn2.*cn1.*dn1)./(cn2.^2 + m*sn1.^2.*sn2.^2);
  yp(:) = log(sn)/pi;
  yprime(:) = sqrt((1-sn.^2)).*sqrt((1-m*sn.^2))./(pi*sn);
else
  high = imag(zp) > Kp/2;
  yp(~high) = (-i*zp(~high) + log(-i/2*(exp(2*i*zp(~high))-1)))/pi;
  yprime(~high) = i*(2./(1-exp(-2*i*zp(~high)))-1)/pi;
  u = i*Kp-zp(high);
  yp(high) = L + i+ (i*u - log(-i/2*(exp(2*i*u)-1)))/pi;
  yprime(high) = i*(2./(1-exp(-2*i*u))-1)/pi;
end

% Make sure everything is in the strip (roundoff could put it outside)
yp = real(yp) + i*max(0,imag(yp));
yp = real(yp) + i*min(1,imag(yp));

