function [p,id]=OrderPoints(p,cm);
% ORDERPOINTS order points counterclockwise
%   [p,id]=OrderPoints(p,cm); orders the points in a plane in 3d
%   stored columnwise in the matrix p counterclockwise looking from
%   the side opposite to the point cm outside the plane. There must
%   be more than two points. p contains the reorderd points, and id
%   contains the index reordering.

d1=p(:,2)-p(:,1);               % two reference vectors
d2=p(:,3)-p(:,1);
if (cm-p(:,1))'*cross(d1,d2)<0  % good orientation ?
  dr=d1/norm(d1);               % 'x-axis' unit vector
  dn=d2-d2'*dr*dr;
  dn=dn/norm(dn);               % 'y-axis' unit vector
else
  dr=d2/norm(d2);               % 'x-axis' unit vector
  dn=d1-d1'*dr*dr;
  dn=dn/norm(dn);               % 'y-axis' unit vector
end;
for j=2:size(p,2)
  d=p(:,j)-p(:,1);
  ao(j-1)=angle(d'*dr+sqrt(-1)*d'*dn);
end;
[tmp,id]=sort(ao);
id=[1 id+1];
p=p(:,id);
