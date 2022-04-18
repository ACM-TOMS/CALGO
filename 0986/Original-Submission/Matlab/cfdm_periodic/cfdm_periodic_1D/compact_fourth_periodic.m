function[D]=compact_fourth_periodic(n,p,x_l,x_r)
% M. Mehra & K. S. Patel
% Inputs
% p=order of accuracy (p=4 or 6)
% n=number of grid points 
% x_l=left end of the interval (Default value is 0)
% x_r=right end of the interval (Default value is 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
% D = differentiation matrix of order n for fourth derivative approximation (d^4/dx^4) when
% periodic boundary coditions are given
if nargin < 3
    x_l=0;x_r=1;
end;
  h=(x_r-x_l)/n;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if p==4
    d=1/4;a=3/2;
idm=ones(n,1);
H1=spdiags([d*idm 1*idm d*idm], -1:1, n,n);
H2=spdiags([1*idm -4*idm 6*idm -4*idm 1*idm], -2:2, n,n);
%............ Incorporating boundary conditions.....
H1(1,1)=1;H1(1,2)=d;H1(1,n)=d;
H1(n,1)=d;H1(n,n)=1;H1(n,n-1)=d;
H2(1,1)=6;H2(1,2)=-4;H2(1,3)=1;H2(1,n)=-4;H2(1,n-1)=1;
H2(2,1)=-4;H2(2,2)=6;H2(2,3)=-4;H2(2,4)=1;H2(2,n)=1;
H2(n,1)=-4;H2(n,2)=1;H2(n,n)=6;H2(n,n-1)=-4;H2(n,n-2)=1;
H2(n-1,1)=1;H2(n-1,n)=-4;H2(n-1,n-1)=6;H2(n-1,n-2)=-4;H2(n-1,n-3)=1;
H2=(a/(h^4))*H2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if p==6
d=7/26;a=19/13;b=1/13;
a=a/(h^4);b=b/(6*h^4);
idm=ones(n,1);
H1=spdiags([d*idm 1*idm d*idm], -1:1, n,n);
H2=spdiags([b*idm a*idm -(4*a+9*b)*idm (6*a+16*b)*idm -(4*a+9*b)*idm a*idm b*idm], -3:3, n,n);
%............ Incorporating boundary conditions.....
H1(1,1)=1;H1(1,2)=d;H1(1,n)=d;
H1(2,1)=d;H1(2,2)=1;H1(2,3)=d;
H1(n,1)=d;H1(n,n)=1;H1(n,n-1)=d;
H1(n-1,n)=d;H1(n-1,n-1)=1;H1(n-1,n-2)=d;
H2(1,1)=(6*a+16*b);H2(1,2)=-(4*a+9*b);H2(1,3)=a;H2(1,4)=b;H2(1,n)=-(4*a+9*b);H2(1,n-1)=a;H2(1,n-2)=b;
H2(2,1)=-(4*a+9*b);H2(2,2)=(6*a+16*b);H2(2,3)=-(4*a+9*b);H2(2,4)=a;H2(2,5)=b;H2(2,n)=a;H2(2,n-1)=b;
H2(3,1)=a;H2(3,2)=-(4*a+9*b);H2(3,3)=(6*a+16*b);H2(3,4)=-(4*a+9*b);H2(3,5)=a;H2(3,6)=b;H2(3,n)=b;
H2(n,1)=-(4*a+9*b);H2(n,2)=a;H2(n,3)=b;H2(n,n)=(6*a+16*b);H2(n,n-1)=-(4*a+9*b);H2(n,n-2)=a;H2(n,n-3)=b;
H2(n-1,1)=a;H2(n-1,2)=b;H2(n-1,n)=-(4*a+9*b);H2(n-1,n-1)=(6*a+16*b);H2(n-1,n-2)=-(4*a+9*b);H2(n-1,n-3)=a;H2(n-1,n-4)=b;
H2(n-2,1)=b;H2(n-2,n)=a;H2(n-2,n-1)=-(4*a+9*b);H2(n-2,n-2)=(6*a+16*b);H2(n-2,n-3)=-(4*a+9*b);H2(n-2,n-4)=a;H2(n-2,n-5)=b;
end
D=H1\H2;
end