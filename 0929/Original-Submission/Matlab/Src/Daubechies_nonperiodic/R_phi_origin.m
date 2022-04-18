function [y1,ro]=R_phi_origin(k,D,h,supp,phi)
%[y1,ro]=R_phi_origin(k,D,h,supp,phi)
% To find the value of left hand boundary function at 0
% Also calulates L_ro(k,k)
% Dependencies
% R_daubfilt.m, R_phi.m, 

M=D/2;   % Here M is number of vanishing moment
[hk_R]=R_daubfilt(k,D);
c=(1-sqrt(2)*hk_R(1))/(sqrt(2)*hk_R(2));
if(D==4)
x(1)=-1/2;
x(2)=-1;
y(1)=R_phi(k,x(1),D,h,supp,phi);
y(2)=R_phi(k,x(2),D,h,supp,phi);
p=polyfit(x,y,M-1);
x1=0;
y1=polyval(p,x1);
%y=y1+((y2-y1)/(x2-x1))*(x-x1); % To find the value of boundary function at 0
elseif(D==8)
    x(1)=-1/8;
    x(2)=-1/4;
    x(3)=-1/2;
    x(4)=-1;
    y(1)=R_phi(k,x(1),D,h,supp,phi);
    y(2)=R_phi(k,x(2),D,h,supp,phi);
    y(3)=R_phi(k,x(3),D,h,supp,phi);
    y(4)=R_phi(k,x(4),D,h,supp,phi);
    p=polyfit(x,y,M-1);
    x1=0;
    y1=polyval(p,x1);
end
ro=(y1^2)/2;  %sign is not there due to definetion of right hand side
