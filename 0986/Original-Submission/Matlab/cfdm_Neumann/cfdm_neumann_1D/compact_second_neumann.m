function[D,K]=compact_second_neumann(n,u_left,u_right,x_l,x_r)
% M. Mehra & K. S. Patel
% Inputs
% n=number of grid points 
% u_left=u_x(0) i.e. left Neumann boundary condition
% u_right=u_x(1)  i.e. right Neumann boundary condition
% x_l=left end of the interval (Default value is 0)
% x_r=right end of the interval (Default value is 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remark: This approximation is fourth order accurate 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
% D= differentiation matrix of order n for second derivative approximation (d^2/dx^2) when Neumann boundary coditions are given
% K= It is a vector of size n x 1 consiting the boundary conditions % f''=D2*f+K2.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 5
    x_l=0;x_r=1;
 end;
h=(x_r-x_l)/(n-1);
idm=ones(n,1);
H1=spdiags([1*idm 10*idm 1*idm], -1:1, n,n);
H2=spdiags([-6*idm 12*idm -6*idm], -1:1, n,n);
H1(1,1)=22;H1(1,2)=-4;H1(n,n)=22;H1(n,n-1)=-4;H2(1,1)=6;H2(1,2)=-6;H2(n,n)=6;H2(n,n-1)=-6;
S1(1,1)=-u_left/h;S1(n,1)=u_right/h;S1=sparse(S1);
A=H1\H2;
D=(-2/(h^2))*A;
B=H1\S1;
K=12*B;
end