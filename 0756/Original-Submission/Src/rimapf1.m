function zdot = rimapf1(wp,yp);
%RIMAPF1 (not intended for calling directly by the user)
%	Used by RINVMAP for solution of an ODE.

global SCIMDATA 

lenyp = length(yp);
lenzp = lenyp/2;
zp = yp(1:lenzp)+sqrt(-1)*yp(lenzp+1:lenyp);
n = SCIMDATA(1,5);

f = SCIMDATA(1:lenzp,1)./rderiv(zp,SCIMDATA(1:n,2),SCIMDATA(1:n,3),...
    SCIMDATA(2,5),SCIMDATA(1:n,4));
zdot = [real(f);imag(f)];
