function[D]=second_compact_dirichlet(n,p,x_l,x_r)
% M. Mehra & K. S. Patel
% Inputs
% n=number of grid points 
% p=order of accuracy (p=4 or 6)
% x_l=left end of the interval (Default value is 0)
% x_r=right end of the interval (Default value is 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
% D= differentiation matrix of order n for second derivative approximation (d^2/dx^2) when
% Dirichlet boundary ccoditions are given
 if nargin < 3
    x_l=0;x_r=1;
  end;
h=(x_r-x_l)/(n-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[a,b,c,d,e]=coefficient_periodic_second(p); % Coefficients for Interior points
a=a/(h*h);b=b/(4*h*h);c=c/(9*h*h);
idm=ones(n,1);
H1=spdiags([d*idm 1*idm d*idm], -1:1, n,n);
[a1,a2,a3,a4,a5,a6,a7,a8,a9,alp]=coefficient_dirichlet_second1(p); % Coefficient for i=1 and N
if p==4
H2=spdiags([a*idm -2*a*idm a*idm], -1:1, n,n);
%.........Incorporating boundary conditions...................
H1(1,1)=1;H1(1,2)=alp;H1(n,n)=1;H1(n,n-1)=alp;
H2(1,1)=a1/(h*h);H2(1,2)=a2/(h*h);H2(1,3)=a3/(h*h);H2(1,4)=a4/(h*h);H2(1,5)=a5/(h*h);
H2(n,n)=a1/(h*h);H2(n,n-1)=a2/(h*h);H2(n,n-2)=a3/(h*h);H2(n,n-3)=a4/(h*h);H2(n,n-4)=a5/(h*h);
%.............................................................
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if p==6
H2=spdiags([b*idm a*idm -2*(a+b)*idm a*idm b*idm], -2:2, n,n);
[b1,b2,b3,b4,b5,b6,b7,b8,b9,bet]=coefficient_dirichlet_second2(p); % Coefficient for i=2 and N-1
%....................Incorporating boundary conditions......
H1(1,1)=1;H1(1,2)=alp;H1(n,n)=1;H1(n,n-1)=alp;
H1(2,1)=bet;H1(2,2)=1;H1(2,3)=bet;H1(n-1,n)=bet;H1(n-1,n-1)=1;H1(n-1,n-2)=bet;
H2(1,1)=a1/(h*h);H2(1,2)=a2/(h*h);H2(1,3)=a3/(h*h);H2(1,4)=a4/(h*h);H2(1,5)=a5/(h*h);H2(1,6)=a6/(h*h);H2(1,7)=a7/(h*h);
H2(n,n)=a1/(h*h);H2(n,n-1)=a2/(h*h);H2(n,n-2)=a3/(h*h);H2(n,n-3)=a4/(h*h);H2(n,n-4)=a5/(h*h);H2(n,n-5)=a6/(h*h);H2(n,n-6)=a7/(h*h);
H2(2,1)=b1/(h*h);H2(2,2)=b2/(h*h);H2(2,3)=b3/(h*h);H2(2,4)=b4/(h*h);H2(2,5)=b5/(h*h);H2(2,6)=b6/(h*h);H2(2,7)=b7/(h*h);
H2(n-1,n)=b1/(h*h);H2(n-1,n-1)=b2/(h*h);H2(n-1,n-2)=b3/(h*h);H2(n-1,n-3)=b4/(h*h);H2(n-1,n-4)=b5/(h*h);H2(n-1,n-5)=b6/(h*h);H2(n-1,n-6)=b7/(h*h);
%..............................................................................
end
D=H1\H2;
end
