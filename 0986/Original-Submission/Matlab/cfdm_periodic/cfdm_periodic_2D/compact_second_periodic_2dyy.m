function[D]=compact_second_periodic_2dyy(n,p,y_l,y_r)
% M. Mehra & K. S. Patel
% Inputs
% p=order of accuracy (p=4 or 6 or 8 or 10)
% n=number of grid points in y direction 
% y_l=left end of the interval in y direction (Default value is 0)
% y_r=right end of the interval in y direction (Default value is 1)
% Output
% D= differentiation matrix of order n for second partial derivative approximation with respect to y (del^2/dely^2) when
% periodic boundary coditions are given
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 if nargin < 3
    y_l=0;y_r=1;
  end;  
nn=n+1;
 [D1]=compact_second_periodic(n,p,y_l,y_r);
%.. Insert this matrix on the diagonal................
 for i=1:n
 D11((nn-1)*(i-1)+1:(nn-1)*i,(nn-1)*(i-1)+1:(nn-1)*i)=D1;
 end
%.........Change this matrix in y direction...
m=1;
for k=1:n
   g=1; 
for i=1:n
    for j=1:n
        D(m,((nn-1)*(j-1))+g)=D11(k,j);
    end
    m=m+1;
   g=g+1;
end
end
%...............................................
end