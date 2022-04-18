function PlotMesh(N,T,col);
% PLOTMESH plots a triangluar mesh
%   PlotMesh(N,T); plots the mesh given by the nodes N and triangles
%   T. The real boundaries are drawn in bold and for small meshes
%   the node numbers are added as well. 

axis('equal');
for i=1:size(T,1),
  for j=1:3,
    line([N(1,T(i,j)) N(1,T(i,mod(j,3)+1))], ...
      [N(2,T(i,j)) N(2,T(i,mod(j,3)+1))],'Color',col);
    bc=mean(N(:,T(i,1:3))')';
  end;
end;
m=size(N,2);
if m<10,
  for i=1:m,
    text(N(1,i)+.01,N(2,i)+.02,num2str(i));
  end;
end;
