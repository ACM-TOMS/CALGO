function zdot = deimapfun(wp,yp,flag,scale,z,beta,c);
%   Used by DEINVMAP for solution of an ODE.

%   Copyright 1998 by Toby Driscoll.
%   $Id: deimapfun.m,v 2.1 1998/05/10 04:33:25 tad Exp $

lenyp = length(yp);
lenzp = lenyp/2;
zp = yp(1:lenzp)+sqrt(-1)*yp(lenzp+1:lenyp);

f = scale./dederiv(zp,z,beta,c);
zdot = [real(f);imag(f)];
