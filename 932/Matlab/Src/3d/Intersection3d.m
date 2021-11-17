function [P,nc,H,v]=Intersection3d(X,Y);
% INTERSECTION3D computes intersection polyhedron of two tetrahedra
%   [P,nc,H,v]=Intersection3d(X,Y); computes for the two given tetrahedra
%   X and Y (point coordinates are stored column-wise) the points P
%   where they intersect, in nc the bolean information of which
%   neighbors of X are also intersecting with Y, in H stored row-wise
%   the polygons of the faces of the intersection polyhedron, stored
%   counter-clock-wise looking from the outside, and in v the volume
%   of the intersection.

P=[];nc=[0 0 0 0];
sxk=zeros(4,1);                        % intersection node index
syk=zeros(4,1);                        % on surfaces of X and Y
l1=[1 2 3 4 1 2];                      % enumeration of lines
l2=[2 3 4 1 3 4];
s1=[1 1 2 3 1 2];                      % faces touching each line
s2=[4 2 3 4 3 4];
ni=[1 3 4; 1 2 4; 1 2 3; 2 3 4];       % faces touching each node
for i=1:6                              % find intersections of edges of
  for j=1:4                            % X with surfaces of Y
    p=TriangleLineIntersection([X(:,l1(i)) X(:,l2(i))],...
      [Y(:,j) Y(:,mod(j,4)+1) Y(:,mod(j+1,4)+1)]);
    if ~isempty(p)
      [k,P]=InsertPoint(p,P);          % insert point if new
      nc(s1(i))=1;nc(s2(i))=1;
      syk(j)=syk(j)+1;SY(j,syk(j))=k;  % remember to which surfaces
      sxk(s1(i))=sxk(s1(i))+1;SX(s1(i),sxk(s1(i)))=k;  % it belongs
      sxk(s2(i))=sxk(s2(i))+1;SX(s2(i),sxk(s2(i)))=k;
    end;
  end;
end;
for i=1:6                              % find intersections of edges of
  for j=1:4                            % Y with surfaces of X
    p=TriangleLineIntersection([Y(:,l1(i)) Y(:,l2(i))],...
      [X(:,j) X(:,mod(j,4)+1) X(:,mod(j+1,4)+1)]);
    if ~isempty(p)
      [k,P]=InsertPoint(p,P);
      nc(j)=1;
      sxk(j)=sxk(j)+1;SX(j,sxk(j))=k;
      syk(s1(i))=syk(s1(i))+1;SY(s1(i),syk(s1(i)))=k;
      syk(s2(i))=syk(s2(i))+1;SY(s2(i),syk(s2(i)))=k;
    end;
  end;
end;

fn=[0 0 0 0];
for i=1:4                              % find interior points of X in Y
  if PointInTetrahedron(X(:,i),Y)
    [k,P]=InsertPoint(X(:,i),P);
    sxk(ni(i,1))=sxk(ni(i,1))+1;SX(ni(i,1),sxk(ni(i,1)))=k;
    sxk(ni(i,2))=sxk(ni(i,2))+1;SX(ni(i,2),sxk(ni(i,2)))=k;
    sxk(ni(i,3))=sxk(ni(i,3))+1;SX(ni(i,3),sxk(ni(i,3)))=k;
    p=X(:,i); fn(i)=1;
  end;
  if PointInTetrahedron(Y(:,i),X)       % and vice versa
    [k,P]=InsertPoint(Y(:,i),P);
    syk(ni(i,1))=syk(ni(i,1))+1;SY(ni(i,1),syk(ni(i,1)))=k;
    syk(ni(i,2))=syk(ni(i,2))+1;SY(ni(i,2),syk(ni(i,2)))=k;
    syk(ni(i,3))=syk(ni(i,3))+1;SY(ni(i,3),syk(ni(i,3)))=k;
    p=Y(:,i);
  end;
end;
if sum(fn)>1                           % more than one point inside
  nc=[1 1 1 1];                        % means all neighbors intersect
end;
H=[];                                  % contains surfaces of the
n=0;                                   % intersection polyhedron
v=0;
if size(P,2)>3                         % construct intersection polyhedra,
  cm=sum(P')'/size(P,2);               % center of mass
  for i=1:4                            % scan surfaces of X and Y
    if sxk(i)>0
      no=RemoveDuplicates(SX(i,1:sxk(i)));
      if length(no)>2 & NewFace(H,no)  % don't compute degenerate polygons
        p=P(:,no);                     % and check if face is new
        [p,id]=OrderPoints(p,cm);      % order points counter-clock-wise
        n=n+1;
        H(n,1:length(id))=no(id);      % add surface
        v=v+PyramidVolume(p,cm);       % add volume
      end;
    end
    if syk(i)>0
      no=RemoveDuplicates(SY(i,1:syk(i)));
      if length(no)>2 & NewFace(H,no)
        p=P(:,no);
        [p,id]=OrderPoints(p,cm);
        n=n+1;
        H(n,1:length(id))=no(id);
        v=v+PyramidVolume(p,cm);
      end;
    end
  end
end;
