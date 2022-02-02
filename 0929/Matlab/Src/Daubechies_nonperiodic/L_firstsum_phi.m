function [L_ph]=L_firstsum_phi(k,x,hk_L,D,M,h,supp,phi)
%[L_ph]=L_firstsum_phi1(k,x,hk_L,D,M,h,supp,phi)
% A function called by L_phi.m and is calling L_phi.m
% Dependencies
% L_phi.m
sum1=0;
for l=0:M-1
     sup=(2*x);
     if(sup>=0 & sup<=M+l)
     sum1=sum1+hk_L(l+1)*L_phi(l,sup,D,h,supp,phi);
     end
end
L_ph=sum1;

