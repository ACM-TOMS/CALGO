function zdot = imapf1(wp,yp);
%HPIMAPF1 (not intended for calling directly by the user)
%	Used by HPINVMAP for solution of an ODE.

global HPIMDATA 

lenyp = length(yp);
lenzp = lenyp/2;
zp = yp(1:lenzp)+sqrt(-1)*yp(lenzp+1:lenyp);
lenx = HPIMDATA(1,4);
bigx = HPIMDATA(1:lenx,2)*ones(1,lenyp/2);
bigbeta = HPIMDATA(1:lenx,3)*ones(1,lenyp/2);

f = HPIMDATA(1:lenzp,1).*exp(sum(log(ones(lenx,1)*zp.' - bigx).*...
    (-bigbeta))).';
zdot = [real(f);imag(f)];
