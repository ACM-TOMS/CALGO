function zdot = dimapf1(wp,yp);
%DIMAPF1 (not intended for calling directly by the user)
%	Used by DINVMAP for solution of an ODE.

global SCIMDATA 

lenyp = length(yp);
lenzp = lenyp/2;
zp = yp(1:lenzp)+sqrt(-1)*yp(lenzp+1:lenyp);
n = SCIMDATA(1,4);
bigz = SCIMDATA(1:n,2)*ones(1,lenzp);
bigbeta = SCIMDATA(1:n,3)*ones(1,lenzp);

f = SCIMDATA(1:lenzp,1).*exp(sum(log(1 - (ones(n,1)*zp.' )./bigz).*...
    (-bigbeta))).';
zdot = [real(f);imag(f)];
