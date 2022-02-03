function [P,n]=EdgeIntersections(X,Y)
% EDGEINTERSECTIONS computes edge intersections of two triangles
%   [P,n]=EdgeIntersections(X,Y) computes for the two given triangles X
%   and Y the points P where their edges intersect. In addition,
%   in n the indices of which neighbors of X are also intersecting
%   with Y are given.

P=[]; k=0; n=[0 0 0];
for i=1:3                            % find all intersections of edges
  for j=1:3
     b=Y(:,j)-X(:,i);
     A=[X(:,mod(i,3)+1)-X(:,i) -Y(:,mod(j,3)+1)+Y(:,j)];
     if rank(A)==2                   % edges not parallel
       r=A\b;
       if r(1)>=0 & r(1)<=1 & r(2)>=0 & r(2)<=1,  % intersection found
         k=k+1; P(:,k)=X(:,i)+r(1)*(X(:,mod(i,3)+1)-X(:,i)); n(i)=1;
       end;
     end;
  end;
end;
