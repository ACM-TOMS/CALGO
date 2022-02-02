function [coef,varargout]=reference(varargin)

% REFERENCE Computes the Argyris basis
%
% COEF=REFERENCE Returns in COEF the coefficients in the mononial 
% basis of the Argyris basis functions in the reference triangle.  
%
% COEF=REFERENCE(P) Returns in COEF the coefficients of the Argyris 
% basis in the triangle specified by P. 
% P must be a 2x3 matrix with P(:,j) the coordinates of vertex j
%
% [COEF,Q]=REFERENCE(P) returns in Q the basis in the triangle 
% specified  by P. Q is a symbolic matrix with Q(i) the ith element 
% of the Argyris basis
%
% NEEDS THE SYMBOLIC TOOLBOX
%

syms x y               % x,y are symbolic

if nargin==0
    p=[0 1 0; 0 0 1];  % reference triangle
else
    p=varargin{1};
    if ~isequal(size(p), [2 3])
        error('First argument must be a 2x3 matrix')
    end
end
coef=[];
coef(1:3,:)=[monomials(p(1,1),p(2,1)); monomials(p(1,2),p(2,2));...
          monomials(p(1,3),p(2,3))];

gradient=[diff(monomials(x,y),x); diff(monomials(x,y),y)];
coef=[coef; subs(gradient,{x,y},{p(1,1),p(2,1)});...
        subs(gradient,{x,y},{p(1,2),p(2,2)});...
        subs(gradient,{x,y},{p(1,3),p(2,3)})];

second=[diff(gradient,x);diff(gradient,y)];
second(3,:)=[]; % ignoring the yx derivative
coef=[coef;subs(second,{x,y},{p(1,1),p(2,1)});...
    subs(second,{x,y},{p(1,2),p(2,2)});...
    subs(second,{x,y},{p(1,3),p(2,3)})];

pm=(p(:,[1 1 2])+p(:,[2 3 3]))/2;  % mid points of each side
sides=(p(:,[2 3 3])-p(:,[1 1 2])); % sides
normal=[-sides(2,:); sides(1,:)];
l=1./sqrt(dot(sides,sides));
normal=normal*diag(l);             % unit normals

normal12=normal(:,1)'*subs(gradient,{x,y},{pm(1,1), pm(2,1)});
normal13=normal(:,2)'*subs(gradient,{x,y},{pm(1,2), pm(2,2)});
normal23=normal(:,3)'*subs(gradient,{x,y},{pm(1,3), pm(2,3)});
coef=[coef;normal12;normal13;normal23];
coef=inv(coef');

% return, if demanded, the basis in symbolic format
if nargout==2
    varargout{1}=coef*monomials(x,y).'; 
end

return

function z=monomials(x,y)
x=x(:); y=y(:);
z=[ x.^0, ...
       x,        y, ...
    x.^2,     x.*y,        y.^2,...
    x.^3,  x.^2.*y,     x.*y.^2,        y.^3,...
    x.^4,  x.^3.*y,  x.^2.*y.^2,     x.*y.^3,     y.^4,...
    x.^5,  x.^4.*y,  x.^3.*y.^2,  x.^2.*y.^3,  x.*y.^4,   y.^5];
