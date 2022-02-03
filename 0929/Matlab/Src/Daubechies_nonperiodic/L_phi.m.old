function [L_phi]=L_phi(k,x,D,h,supp,phi)
% [L_phi]=L_phi(k,x,D,h,supp,phi)
% This function Calculates the value of left boundary scaling function
% Dependencies
%L_daubfilt.m, L_firstsum_phi.m
M=D/2 ;  % Here M is number of vanishing moment
[hk_L]=L_daubfilt(k,D);
     [L_ph]=L_firstsum_phi(k,x,hk_L,D,M,h,supp,phi);
     sum1=L_ph;
sum2=0;
% Support of phi is from [-M+1,M]................
for m=M:M+2*k
    2*x-m;
    ind=(-(-M+1)+(2*x-m))*(1/h)+1;
    if(ind>=1 && ind<=supp)
        hk_L(m+1);
    ph=phi(ind);
    sum2=sum2+hk_L(m+1)*phi(ind);
    end
end
sum=sqrt(2)*(sum1+sum2);
L_phi=sum;
