function [dxx,dxy,dyy]= evalHessArgy(x,y,p)

% EVALHESSARGY returns the second derivatives of the Argyris basis
%
% [DXX,DXY,DYY]=EVALHESSARGY(X,Y,P) returns in DXX, DXY, DYY the 
% values  of the XX-, XY- and YY- derivatives of the Argyris basis on 
% the triangle given by P evaluated at (X,Y)
% 
% P is a 2 x 3 matrix with the coordinates of the vertices of the 
% triangle. X,Y can be row vectors.
%
% DXX,DXY, DYY are 21 x m matrices with 
%
% DXX(i,j)=\partial_{XX} N_i(X(j),Y(j)) 
% DXY(i,j)=\partial_{XY} N_i(X(j),Y(j))
% DYY(i,j)=\partial_{YY} N_i(X(j),Y(j))  i=1...21, j=1...m
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

[C,B,b,Th]=changeOfBasis(p);
[x,y]=k2khat(x,y,B,b);
k=length(x);
mxx=derxx(x,y);
mxy=derxy(x,y);
myy=deryy(x,y);
hessian=zeros(21,3*k);
hessian(:)=[mxx(:) mxy(:) myy(:)]/Th';
hessian=C'*coefRef*hessian;
dxx=hessian(:,1:k);
dxy=hessian(:,k+1:2*k);
dyy=hessian(:,2*k+1:3*k);

return


function z=derxx(x,y)

z=[      0*x;         ...
         0*x;        0*y;         ...
      2*x.^0;        0*y;        0*y;       ...
         6*x;        2*y;        0*y;      0*y;    ...
     12*x.^2;     6*x.*y;     2*y.^2;   0*y.^3;   0*y;  ...
     20*x.^3; 12*x.^2.*y;  6*x.*y.^2;   2*y.^3;   0*y;  0*y];
return


function z=derxy(x,y)

z=[    0*x;        ... 
       0*x;    0*x.^0;         ...
       0*x;      x.^0;        0*x;          ...      
       0*x;       2*x;        2*y;         0*x;      ...
       0*x;    3*x.^2;     4*x.*y;      3*y.^2;     0*x;   ...
       0*x;    4*x.^3;  6*x.^2.*y;   6*x.*y.^2;  4*y.^3;  0*x];
return


function z=deryy(x,y)

z=[    0*x;      ...
       0*x;     0*y;        ...
       0*x;     0*x;    2*y.^0;         ...
       0*x;     0*x;       2*x;        6*y;          ...
       0*x;     0*x;    2*x.^2;     6*x.*y;     12*y.^2;      ...
       0*x;     0*x;    2*x.^3;  6*x.^2.*y;  12*x.*y.^2;  20*y.^3];
return


