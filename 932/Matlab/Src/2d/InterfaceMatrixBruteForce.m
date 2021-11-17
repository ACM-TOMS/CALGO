function M=InterfaceMatrixBruteForce(Na,Ta,Nb,Tb);
% INTERFACEMATRIXBRUTEFORCE projection matrix for nonmatching grids
%   M=InterfaceMatrixBruteForce(Na,Ta,Nb,Tb); takes two triangular
%   meshes Ta and Tb with associated nodal coordinates in Na and Nb
%   and computes the interface projection matrix M using a brute
%   force approach

M=sparse(size(Nb,2),size(Na,2));
for i=1:size(Ta,1)
  for j=1:size(Tb,1)
    [P,nc,Mc]=Intersect(Nb(:,Tb(j,1:3)),Na(:,Ta(i,1:3)));
    if ~isempty(P)             % intersection found
      M(Tb(j,1:3),Ta(i,1:3))=M(Tb(j,1:3),Ta(i,1:3))+Mc;
    end
  end
end
