function[D]=compact_second_periodic_2dxx(n,p,x_l,x_r)
% M. Mehra & K. S. Patel
% Inputs
% n=number of grid points in x direction
% p=order of accuracy (p=4 or 6 or 8 or 10)
% x_l=left end of the interval in x direction (Default value is 0)
% x_r=right end of the interval in x direction (Default value is 1)
% Output
% D = differentiation matrix of order n for second partial derivative approximation with respect to x (del^2/delx^2) when
% periodic boundary coditions are given
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 if nargin < 3
    x_l=0;x_r=1;
  end;  
nn=n+1;
[D1]=compact_second_periodic(n,p,x_l,x_r);
%.. Insert this matrix on the diagonal..........
for i=1:n
D((nn-1)*(i-1)+1:(nn-1)*i,(nn-1)*(i-1)+1:(nn-1)*i)=D1;
end
%...........................................
end
