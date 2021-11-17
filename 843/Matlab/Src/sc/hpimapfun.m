function zdot = hpimapfun(wp,yp,flag,scale,z,beta,c);
%   Used by HPINVMAP for solution of an ODE.

%   Copyright 1998 by Toby Driscoll.
%   $Id: hpimapfun.m,v 2.2 2001/02/21 17:08:53 driscoll Exp $

lenyp = length(yp);
lenzp = lenyp/2;

% Don't allow points in lower half-plane. This really messes up the
% derivative calculation.
zp = yp(1:lenzp) + i*max(0,yp(lenzp+1:lenyp));

f = scale./hpderiv(zp,z,beta,c);
zdot = [real(f);imag(f)];
