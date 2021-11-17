function [y1,ro]=L_phi_origin(k,D,h,supp,phi)
%[y,ro]=L_phi_origin(k,D,h,supp,phi)
% To find the value of left hand boundary function at 0
% Also calulating L_ro(k,k)
%Function called by L_ro.m
% Dependencies
% L_daubfilt.m, L_phi.m


M=D/2;   % Here M is number of vanishing moment
[hk_L]=L_daubfilt(k,D);
if(D==4)
x(1)=1/2;
x(2)=1;
y(1)=L_phi(k,x(1),D,h,supp,phi);
y(2)=L_phi(k,x(2),D,h,supp,phi);
p=polyfit(x,y,M-1);
x1=0;
y1=polyval(p,x1);
elseif(D==8)
    x(1)=1/8;
    x(2)=1/4;
    x(3)=1/2;
    x(4)=1;
    y(1)=L_phi(k,x(1),D,h,supp,phi);
    y(2)=L_phi(k,x(2),D,h,supp,phi);
    y(3)=L_phi(k,x(3),D,h,supp,phi);
    y(4)=L_phi(k,x(4),D,h,supp,phi);
    p=polyfit(x,y,M-1);
    x1=0;
    y1=polyval(p,x1);

end
ro=-(y1^2)/2;
