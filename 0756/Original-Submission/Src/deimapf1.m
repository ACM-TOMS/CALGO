function zdot = deimapf1(wp,yp);
%DEIMAPF1 (not intended for calling directly by the user)
%	Used by DEINVMAP for solution of an ODE.

global SCIMDATA 

lenyp = length(yp);
lenzp = lenyp/2;
zp = yp(1:lenzp)+sqrt(-1)*yp(lenzp+1:lenyp);
lenz = SCIMDATA(1,4);
bigz = SCIMDATA(1:lenz,2)*ones(1,lenzp);
bigbeta = SCIMDATA(1:lenz,3)*ones(1,lenzp);

f = SCIMDATA(1:lenzp,1).*exp(sum(log(1 - (ones(lenz,1)*zp.' )./bigz).*...
    (-bigbeta))).'.*zp.^2;
zdot = [real(f);imag(f)];
