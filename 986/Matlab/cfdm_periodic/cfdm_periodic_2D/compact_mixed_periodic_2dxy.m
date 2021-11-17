function[D]=compact_mixed_periodic_2dxy(n,p,x_l,x_r,y_l,y_r)
% M. Mehra & K. S. Patel
% Inputs
% n=number of grid points in x direction = number of grid points in y direction
% p=order of accuracy (p=4 or 6 or 8 or 10)
% x_l=left end of the interval in x direction (Default value is 0)
% x_r=right end of the interval in x direction (Default value is 1)
% y_l=left end of the interval in y direction (Default value is 0)
% y_r=right end of the interval in y direction (Default value is 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
% D= differentiation matrix of order n for mixed derivative approximation (d^2/dxdy) when
% periodic boundary coditions are given
 if nargin < 3
    x_l=0;x_r=1;y_l=0;y_r=1;
  end;
[E1]=compact_first_periodic(n,p,x_l,x_r);
[E2]=compact_first_periodic(n,p,y_l,y_r);
D=kron(E1,E2);
end