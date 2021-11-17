function[D]=first_compact_dirichlet_2dx(n,p,x_l,x_r)
% M. Mehra & K. S. Patel
% Inputs
% n=number of grid points in x direction
% p=order of accuracy (p=4 or 6)
% x_l=left end of the interval (Default value is 0)
% x_r=right end of the interval (Default value is 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
% D= differentiation matrix of order n for first partial derivative approximation with respect to x (del/delx) when
% Dirichlet boundary ccoditions are given
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 if nargin < 3
    x_l=0;x_r=1;
  end;
[D1]=first_compact_dirichlet(n,p,x_l,x_r);
%.... Insert this matrix on the diagonal.......................
for i=1:n
D((n)*(i-1)+1:(n)*i,(n)*(i-1)+1:(n)*i)=D1;
end
%.............................................................
end
