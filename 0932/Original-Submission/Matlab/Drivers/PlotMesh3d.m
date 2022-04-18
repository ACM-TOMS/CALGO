function PlotMesh3d(N,T,col);
% PLOTMESH3D plots a tetrahedra mesh 
%   PlotMesh3d(N,T,col); plots the mesh given by the nodes N and
%   tetrahedra T in color col. For small meshes the node numbers are
%   added as well.

axis('equal');
view(-37.5+180,30);
for i=1:size(T,1),
  PlotTetrahedron(N(:,T(i,1:4)),col);
end;
m=size(N,2);
de=0.05;
if m<10,
  for i=1:m,
    text(N(1,i)+de,N(2,i)+de,N(3,i)+de,num2str(i),'Color','b');
  end;
end;
