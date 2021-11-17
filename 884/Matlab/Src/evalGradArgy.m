function [dx,dy]=evalGradArgy(x,y,p)

% EVALGRADARGY returns the first derivatives of the Argyris basis
%
% [DX,DY]=EVALGRADARGY(X,Y,P) returns in DX, DY the values of the 
% X- and Y- derivatives of Argyris basis on the triangle given by P 
% evaluated at (X,Y)
% 
% P is a 2 x 3 matrix with the coordinates of the vertices of the 
% triangle. X,Y can be row vectors.
%
% DX,DY are 21 x m matrices with 
%
% DX(i,j)=\partial_X N_i(X(j),Y(j))  
% DY(i,j)=\partial_Y N_i(X(j),Y(j))  i=1...21, j=1..m
%
% m being the length of X, Y
%

global coefRef
if isempty(coefRef) 
    load coefRef.dat    % define coefRef if it does not exist
end
    
if ~isequal(size(p),[2 3])
    error('The third argument must be a 2x3 matrix')
end
dimx=size(x); dimy=size(y);
if ~isequal(dimx,dimy) || length(dimx)>2 || dimx(1)~=1
    error('The first two arguments must be row vectors of the same length')
end

[C,B,b]=changeOfBasis(p);
[x,y]=k2khat(x,y,B,b);
k=length(x);
mx=derx(x,y);
my=dery(x,y);
grads=zeros(21,2*k);
grads(:)=[mx(:) my(:)]/B;
grads=C'*coefRef*grads;
dx=grads(:,1:k);
dy=grads(:,k+1:2*k);

return


function z=derx(x,y)

z=[   0*x; ... 
     x.^0;        0*y; ...
      2*x;          y;           0*y; ...
   3*x.^2;     2*x.*y;          y.^2;        0*y; ...
   4*x.^3;  3*x.^2.*y;     2*x.*y.^2;       y.^3;   0*y; ...
   5*x.^4;  4*x.^3.*y;  3*x.^2.*y.^2;  2*x.*y.^3;  y.^4;   0*y];
return


function z=dery(x,y)

z=[   0*x; ... 
      0*x;   y.^0; ...
      0*x;      x;       2*y; ...
      0*x;   x.^2;    2*x.*y;        3*y.^2; ...
      0*x;   x.^3; 2*x.^2.*y;     3*x.*y.^2;     4*y.^3; ...
      0*x;   x.^4; 2*x.^3.*y;  3*x.^2.*y.^2;  4*x.*y.^3;  5*y.^4];
return