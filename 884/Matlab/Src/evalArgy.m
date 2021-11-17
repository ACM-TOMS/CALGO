function z=evalArgy(x,y,p)

% EVALARGY(X,Y,P) evaluates the Argyris basis
%
% Z=EVALARGY(X,Y,P) returns in Z the values of the Argyris basis 
% on the triangle given by P evaluated at (X,Y)
% 
% P is a 2 x 3 matrix with the coordinates of the vertices of the 
% triangle. X,Y can be row vectors.
% 
% Z is 21 x m matrix with
%
% Z(i,j)=N_i(X(j),Y(j))    i=1...21, m=1..21
%
% m being the length of X, Y 


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
z=monomials(x,y);
z=C'*coefRef*z;

return

function z=monomials(x,y)

z=[x.^0; ... 
      x;        y; ... 
   x.^2;     x.*y;        y.^2; ...
   x.^3;  x.^2.*y;     x.*y.^2;        y.^3; ...
   x.^4;  x.^3.*y;  x.^2.*y.^2;     x.*y.^3;     y.^4; ...
   x.^5;  x.^4.*y;  x.^3.*y.^2;  x.^2.*y.^3;  x.*y.^4;   y.^5];
