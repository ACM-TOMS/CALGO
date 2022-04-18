function [R_phi]=R_phi(k,x,D,h,supp,phi)
% [R_phi]=R_phi(k,x,D,h,supp,phi)
% Calculates the values of right boundary scaling functions
% Dependencies
% R_daubfilt.m, R_firstsum_phi.m, 
M=D/2;   % Here M is number of vanishing moment
[hk_R]=R_daubfilt(k,D);
     [R_ph]=R_firstsum_phi(k,x,hk_R,D,M,h,supp,phi);
     sum1=R_ph;
sum2=0;
for m=M:M+2*k
    2*x+m+1;
   -(-M+1)+(2*x+m+1);  % Here is change from Left hand boundary function.
    ind=(-(-M+1)+(2*x+m+1))*(1/h)+1;
    if(ind>=1 && ind<=supp)
    ph=phi(ind);
    sum2=sum2+hk_R(m+1)*phi(ind);
    end
end
sum=sqrt(2)*(sum1+sum2);
R_phi=sum;
