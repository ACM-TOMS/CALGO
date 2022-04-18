function P=PointsOfXInY(X,Y);
% POINTSOFXINY finds corners of one triangle within another one
%   P=PointsOfXInY(X,Y); computes for the two given triangles X
%   and Y (point coordinates are stored column-wise, in counter clock
%   order) the corners P of X which lie in the interior of Y.

k=0;P=[];
v0=Y(:,2)-Y(:,1); v1=Y(:,3)-Y(:,1);  % find interior points of X in Y
d00=v0'*v0; d01=v0'*v1; d11=v1'*v1;  % using baricentric coordinates
id=1/(d00*d11-d01*d01);
for i=1:3
  v2=X(:,i)-Y(:,1); d02=v0'*v2; d12=v1'*v2;
  u=(d11*d02-d01*d12)*id; v=(d00*d12-d01*d02)*id;
  if u>=0 & v>=0 & u+v<=1            % also include nodes on the boundary
    k=k+1; P(:,k)=X(:,i);
  end;
end;
