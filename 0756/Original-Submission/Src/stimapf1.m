function zdot = stimapf1(wp,yp);
%STIMAPF1 (not intended for calling directly by the user)
%	Used by STINVMAP for solution of an ODE.

global SCIMDATA 

lenyp = length(yp);
lenzp = lenyp/2;
zp = yp(1:lenzp)+sqrt(-1)*yp(lenzp+1:lenyp);
n = SCIMDATA(1,4);

f = SCIMDATA(1:lenzp,1)./stderiv(zp,SCIMDATA(1:n,2),SCIMDATA(1:n,3));
zdot = [real(f);imag(f)];
