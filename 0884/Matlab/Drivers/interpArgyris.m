function [t,s]=interpArgyris(t0,x,y,g,varargin)

%  INTERPARGYRIS Compute the Argyris interporlant 
%
%  [T,S]= INTERPARGYRIS(T0,X0,Y0,G) 
%  computes the Argyris interpolant on the grid specified 
%  by T0,X0,Y0. 
%
%  T0,X0,YO is a triangular grid
%
%  G is an inline function. The derivatives up to order 2 
%  are computed using the symbolic toolbox.
%
%  T is nT x 21 matrices storing the information for the Argyris 
%  elements. Specifically T(i,:) has the global coordinates of 
%  the local degrees of freedom for the Argyris element on triangle i. 
%
%  S is a vector with the coordinates of the Argyris interpolant
%  in the global basis.
%  
%  [T,S]= INTERPARGYRIS(T0,X0,Y0,G,GX,GY,GXX,GXY,GYY) 
%  GX,GY, GXX, GXY, GYY are the derivatives up to order 2 of G. 
%
%  USE IF THE SYMBOLIC TOOLBOX IS NOT INSTALLED
%

f=vectorize(g);             % f is g vectorized
if nargin==4
    % Perform the derivatives
    % Symbolic toolbox is required!!
    syms u v  
    fx=vectorize(inline([char(diff(f(u,v),u)) '.*u.^0'],'u','v'));
    fy=vectorize(inline([char(diff(f(u,v),v)) '.*u.^0'],'u','v'));
    fxx=vectorize(inline([char(diff(fx(u,v),u)) '.*u.^0'],'u','v'));
    fxy=vectorize(inline([char(diff(fy(u,v),u)) '.*u.^0'],'u','v'));
    fyy=vectorize(inline([char(diff(fy(u,v),v)) '.*u.^0'],'u','v'));
elseif nargin==9
    fx=vectorize(varargin{1});
    fy=vectorize(varargin{2});
    fxx=vectorize(varargin{3});
    fxy=vectorize(varargin{4});
    fyy=vectorize(varargin{5});
else
    error('incorrect number of input arguments')
end
    
% Construct the grid where the Argyris function
% is evaluated 
nTri=size(t0,1); 
nPt=length(x);
% Adjacent matrix
adj=sparse(nPt,nPt);
nodos=zeros(nPt,5);
nEl=nPt; % number of element

% add 18 degrees of freedom to each triangle
% we will assign values below
t=[t0 zeros(nTri,18)];

for j=1:nTri
   np=t(j,:);
   taux=zeros(1,18);
   for i=1:3
      if sum(nodos(np(i),:))==0
            nodos(np(i),1:5)=(nEl+1):(nEl+5);
            nEl=nEl+5;
      end
      taux(5*i-4:5*i)=nodos(np(i),:);
   end
   v=np([1 1 2]);
   v2=np([2 3 3]);
   for i=1:3
      if adj(v(i),v2(i))==0
           adj(v(i),v2(i))=(nEl+1);
           adj(v2(i),v(i))=(nEl+1);
            nEl=nEl+1;
      end
      taux(15+i)=abs(adj(v(i),v2(i)));
   end
   % Reordering the nodes 
   % first the x,y derivatives
   % second the xx, xy and yy derivatives
   % third, the normal derivatives at the midpoints of the edges
   aux=taux([1 2 6 7 11 12 3 4 5 8 9 10 13 14 15 16 17 18]);
   t(j,4:21)=aux;   
end


% evaluations of the function and its derivatives up to order 2

eval  = f(x,y);    
difx  = fx(x,y);
dify  = fy(x,y);
difxx = fxx(x,y);
difxy = fxy(x,y);
difyy = fyy(x,y);

s=zeros(1,max(max(t)));
for j=1:size(t,1)
   % values of f at the vertices
   c=t(j,:);
   s(c(1:3))=eval(c(1:3));
   % values of the derivatives of f at the vertices
   for k=1:3
      s(c([4 6 8]))    = difx(c(1:3)); 
      s(c([5 7 9]))    = dify(c(1:3));
      s(c([10 13 16])) = difxx(c(1:3));
      s(c([11 14 17])) = difxy(c(1:3));
      s(c([12 15 18])) = difyy(c(1:3));
   end
   % compute the normal derivatives at the midpoints of the edges
   
   % mid points
   p1x     = (x(c(1))+x(c(2)))/2; p1y=(y(c(1))+y(c(2)))/2;
   p2x     = (x(c(1))+x(c(3)))/2; p2y=(y(c(1))+y(c(3)))/2;
   p3x     = (x(c(2))+x(c(3)))/2; p3y=(y(c(2))+y(c(3)))/2;
   
   % normal vector to the sides (positively rotated about pi/2 )
   w1      = [-(y(c(2))- y(c(1))) x(c(2))- x(c(1))];
   w2      = [-(y(c(3))- y(c(1))) x(c(3))- x(c(1))];
   w3      = [-(y(c(3))- y(c(2))) x(c(3))- x(c(2))];
   
   % Taking into account the orientation of each side...
   %  1 if possitive, 
   % -1 if negative
   m=[c(1)<c(2) c(1)< c(3) c(2)<c(3)];
   m=2*m-1;
   
   % store the values of the normal derivative at the mid points
   % of each side. 
   s(c(19))= (fx(p1x,p1y)*w1(1)+fy(p1x,p1y)*w1(2))/norm(w1)*m(1);
   s(c(20))= (fx(p2x,p2y)*w2(1)+fy(p2x,p2y)*w2(2))/norm(w2)*m(2);
   s(c(21))= (fx(p3x,p3y)*w3(1)+fy(p3x,p3y)*w3(2))/norm(w3)*m(3);
end
return
