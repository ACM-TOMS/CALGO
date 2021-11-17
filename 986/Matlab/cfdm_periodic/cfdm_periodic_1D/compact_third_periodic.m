function[D]=compact_third_periodic(n,p,x_l,x_r)
% M. Mehra & K. S. Patel
% Inputs
% n=number of grid points 
% p=order of accuracy (p =6) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remark: For p=4 results are not accurate since tridiagonal scheme
% becomes singular near pi. It is explained in the manuscript. So we take p=6 only.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x_l=left end of the interval (Default value is 0)
% x_r=right end of the interval (Default value is 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
% D = differentiation matrix of order n for third derivative approximation (d^3/dx^3) when
% periodic boundary coditions are given
 if nargin < 3
    x_l=0;x_r=1;
  end;
  h=(x_r-x_l)/n;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 if p==4
    d=1/2;a=2;a=a/(2*h^3);idm=ones(n,1);
H1=spdiags([d*idm 1*idm d*idm], -1:1, n,n);
H2=spdiags([-a*idm 2*a*idm 0*idm -2*a*idm a*idm], -2:2, n,n);
%............ Incorporating boundary conditions.....
H1(1,1)=1;H1(1,2)=d;H1(1,n)=d;
H1(n,1)=d;H1(n,n)=1;H1(n,n-1)=d;
H2(1,1)=0;H2(1,2)=-2*a;H2(1,3)=a;H2(1,n)=2*a;H2(1,n-1)=-a;
H2(2,1)=2*a;H2(2,2)=0;H2(2,3)=-2*a;H2(2,4)=a;H2(2,n)=-a;
H2(n,1)=-2*a;H2(n,2)=a;H2(n,n)=0;H2(n,n-1)=2*a;H2(n,n-2)=-a;
H2(n-1,1)=a;H2(n-1,n)=-2*a;H2(n-1,n-1)=0;H2(n-1,n-2)=2*a;H2(n-1,n-3)=-a;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if p==6
d=7/16;a=2;b=-1/8;a=a/(2*h^3);b=b/(8*h^3);idm=ones(n,1);
H1=spdiags([d*idm 1*idm d*idm], -1:1, n,n);
H2=spdiags([-b*idm -a*idm (2*a+3*b)*idm 0*idm -(2*a+3*b)*idm a*idm b*idm], -3:3, n,n);
%............ Incorporating boundary conditions.....
H1(1,1)=1;H1(1,2)=d;H1(1,n)=d;
H1(2,1)=d;H1(2,2)=1;H1(2,3)=d;
H1(n,1)=d;H1(n,n)=1;H1(n,n-1)=d;
H1(n-1,n)=d;H1(n-1,n-1)=1;H1(n-1,n-2)=d;
H2(1,1)=0;H2(1,2)=-(2*a+3*b);H2(1,3)=a;H2(1,4)=b;H2(1,n)=(2*a+3*b);H2(1,n-1)=-a;H2(1,n-2)=-b;
H2(2,1)=(2*a+3*b);H2(2,2)=0;H2(2,3)=-(2*a+3*b);H2(2,4)=a;H2(2,5)=b;H2(2,n)=-a;H2(2,n-1)=-b;
H2(3,1)=-a;H2(3,2)=(2*a+3*b);H2(3,3)=0;H2(3,4)=-(2*a+3*b);H2(3,5)=a;H2(3,6)=b;H2(3,n)=-b;
H2(n,1)=-(2*a+3*b);H2(n,2)=a;H2(n,3)=b;H2(n,n)=0;H2(n,n-1)=(2*a+3*b);H2(n,n-2)=-a;H2(n,n-3)=-b;
H2(n-1,1)=a;H2(n-1,2)=b;H2(n-1,n)=-(2*a+3*b);H2(n-1,n-1)=0;H2(n-1,n-2)=(2*a+3*b);H2(n-1,n-3)=-a;H2(n-1,n-4)=-b;
H2(n-2,1)=b;H2(n-2,n)=a;H2(n-2,n-1)=-(2*a+3*b);H2(n-2,n-2)=0;H2(n-2,n-3)=(2*a+3*b);H2(n-2,n-4)=-a;H2(n-2,n-5)=-b;
end
D=H1\H2;
end