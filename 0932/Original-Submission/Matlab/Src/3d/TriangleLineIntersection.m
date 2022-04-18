function p=TriangleLineIntersection(X,Y);
% TRIANGLELINEINTERSECTION intersection of a line and a triangle
%   p=TriangeLineIntersection(X,Y); computes for a given line X in 3d,
%   and a triangle Y in 3d, the point p of intersection in the
%   triangle, and otherwise returns p=[];

p=[];
b=Y(:,1)-X(:,1);
A=[X(:,2)-X(:,1) Y(:,1)-Y(:,2) Y(:,1)-Y(:,3)];
if rank(A)==3                         % edge and triangle not parallel
  r=A\b;
  if r(1)>=0 & r(1)<=1 & r(2)>=0 & r(3)>=0 & r(2)+r(3)<=1
    p=X(:,1)+r(1)*(X(:,2)-X(:,1));
  end;
end;
