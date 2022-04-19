function [P,n,M]=Intersect(X,Y);
% INTERSECT intersection of two triangles and mortar contribution
%   [P,n,M]=Intersect(X,Y); computes for the two given triangles X
%   and Y the points P where they intersect, in n the indices of
%   which neighbors of X are also intersecting with Y, and the
%   local mortar matrix M of contributions of the element X on the 
%   element Y. 

[P,n]=EdgeIntersections(X,Y);
Q=PointsOfXInY(X,Y);
if size(Q,2)>1                       % if two or more interior
  n=[1 1 1];                         % points the triangle is 
end                                  % candidate for all neighbors
P=[P Q];
P=[P PointsOfXInY(Y,X)];
P=SortAndRemoveDoubles(P);           % sort counter clock wise
M=zeros(3,3);
if size(P,2)>0
  for j=2:size(P,2)-1                % compute interface matrix
    M=M+MortarInt(P(:,[1 j j+1]),X,Y);
  end;
%  patch(P(1,:),P(2,:),'m')           % draw intersection for illustration
%  pause(1)
end;
