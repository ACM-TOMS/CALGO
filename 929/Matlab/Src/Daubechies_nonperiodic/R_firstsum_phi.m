function [R_ph]=R_phi1(k,x,hk_R,D,M,h,supp,phi)
% [R_ph]=R_firstsum_phi(k,x,hk_R,D,M,h,supp,phi)
%function called by R_phi.m and is calling R_phi.m
% Dependencies
% R_phi.m
sum1=0;
for l=0:M-1
     sup=(2*x);
     if(sup>=-(M+l) & sup<=0)
     sum1=sum1+hk_R(l+1)*R_phi(l,sup,D,h,supp,phi);
     end
end
R_ph=sum1;

