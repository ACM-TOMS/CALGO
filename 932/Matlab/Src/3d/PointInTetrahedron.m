function p=PointInTetrahedron(X,Y);
% POINTINTETRAHEDRON check if a point is inside a tetrahedron
%   p=PointInTetrahedron(X,Y); checks if the point X is contained in
%   the tetrahedron Y.

D0=det([Y; ones(1,4)]);
D1=det([X Y(:,2:4); ones(1,4)]);
D2=det([Y(:,1) X Y(:,3:4); ones(1,4)]);
D3=det([Y(:,1:2) X Y(:,4); ones(1,4)]);
D4=det([Y(:,1:3) X; ones(1,4)]);
p=sign(D0)==sign(D1) & sign(D0)==sign(D2) & sign(D0)==sign(D3) & ...
  sign(D0)==sign(D4);
