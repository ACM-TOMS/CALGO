function[D]=first_compact_dirichlet_2dy(n,p,y_l,y_r)
% M. Mehra & K. S. Patel
% Inputs
% n=number of grid points in y direction
% p=order of accuracy (p=4 or 6)
% y_l=left end of the interval in y direction (Default value is 0)
% y_r=right end of the interval in y direction (Default value is 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
% D= differentiation matrix of order n for first partial derivative approximation with respect to y (del/dely) when
% Dirichlet boundary ccoditions are given
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 3
    y_l=0;y_r=1;
end;
 [D1]=first_compact_dirichlet(n,p,y_l,y_r);
 %... Insert this matrix on the diagonal......................
for i=1:n
 D11((n)*(i-1)+1:(n)*i,(n)*(i-1)+1:(n)*i)=D1;
end
%......Change this matrix in y-direction.................. 
m=1;
for k=1:n
   g=1; 
for i=1:n
    for j=1:n
        D(m,((n)*(j-1))+g)=D11(k,j);
    end
    m=m+1;
   g=g+1;
end
end
%..........................................................
end