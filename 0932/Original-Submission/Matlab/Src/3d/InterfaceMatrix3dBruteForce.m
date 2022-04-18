function InterfaceMatrix3dBruteForce(Na,Ta,Nb,Tb);
% INTERFACEMATRIX3D projection matrix for nonmatching thetrahedra grids 
%   M=InterfaceMatrix3d(Na,Ta,Nb,Tb); takes two tetrahedra meshes Ta
%   and Tb with associated nodal coordinates in Na and Nb and
%   computes the projection matrix M. As a precondition, the first
%   thetrahedra in Ta and Tb need to intersect.

M=sparse(size(Nb,2),size(Na,2));
for i=1:size(Ta,1)
  for j=1:size(Tb,1)
    [P,nc,H,v]=Intersection3d(Nb(:,Tb(j,1:4)),Na(:,Ta(i,1:4)));
    if ~isempty(P)             % intersection found
    end
  end
end
