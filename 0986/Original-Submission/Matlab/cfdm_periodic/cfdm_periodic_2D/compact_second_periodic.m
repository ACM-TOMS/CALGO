function[D]=compact_second_periodic(n,p,x_l,x_r)
% M. Mehra & K. S. Patel
% Inputs
% n=number of grid points 
% p=order of accuracy (p=4 or 6 or 8 or 10)
% x_l=left end of the interval (Default value is 0)
% x_r=right end of the interval (Default value is 1)
% Output
% D = differentiation matrix of order n for second derivative approximation (d^2/dx^2) when periodic boundary coditions are given
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
if nargin < 3
    x_l=0;x_r=1;
end;
  h=(x_r-x_l)/n;
[a,b,c,d,e]=coefficient_periodic_second(p);
a=a/(h*h);b=b/(4*h*h);c=c/(9*h*h);
idm=ones(n,1);
H1=spdiags([e*idm d*idm 1*idm d*idm e*idm], -2:2, n,n);
H2=spdiags([c*idm b*idm a*idm -(2*a+2*b+2*c)*idm a*idm b*idm c*idm], -3:3, n,n);
%........Incorporating boundary conditions........................
H1(1,1)=1;H1(1,2)=d;H1(1,3)=e;H1(1,n)=d;H1(1,n-1)=e;
H1(2,1)=d;H1(2,2)=1;H1(2,3)=d;H1(2,4)=e;H1(2,n)=e;
H1(n,1)=d;H1(n,2)=e;H1(n,n)=1;H1(n,n-1)=d;H1(n,n-2)=e;
H1(n-1,1)=e;H1(n-1,n)=d;H1(n-1,n-1)=1;H1(n-1,n-2)=d;H1(n-1,n-3)=e;
H2(1,1)=-(2*a+2*b+2*c);H2(1,2)=a;H2(1,3)=b;H2(1,4)=c;H2(1,n)=a;H2(1,n-1)=b;H2(1,n-2)=c;
H2(2,1)=a;H2(2,2)=-(2*a+2*b+2*c);H2(2,3)=a;H2(2,4)=b;H2(2,5)=c;H2(2,n)=b;H2(2,n-1)=c;
H2(3,1)=b;H2(3,2)=a;H2(3,3)=-(2*a+2*b+2*c);H2(3,4)=a;H2(3,5)=b;H2(3,6)=c;H2(3,n)=c;
H2(n,1)=a;H2(n,2)=b;H2(n,3)=c;H2(n,n)=-(2*a+2*b+2*c);H2(n,n-1)=a;H2(n,n-2)=b;H2(n,n-3)=c;
H2(n-1,1)=b;H2(n-1,2)=c;H2(n-1,n)=a;H2(n-1,n-1)=-(2*a+2*b+2*c);H2(n-1,n-2)=a;H2(n-1,n-3)=b;H2(n-1,n-4)=c;
H2(n-2,1)=c;H2(n-2,n)=b;H2(n-2,n-1)=a;H2(n-2,n-2)=-(2*a+2*b+2*c);H2(n-2,n-3)=a;H2(n-2,n-4)=b;H2(n-2,n-5)=c;
D=H1\H2;
end