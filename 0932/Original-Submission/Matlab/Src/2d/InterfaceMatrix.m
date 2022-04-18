function M=InterfaceMatrix(Na,Ta,Nb,Tb);
% INTERFACEMATRIX projection matrix for nonmatching triangular grids
%   M=InterfaceMatrix(Na,Ta,Nb,Tb); takes two triangular meshes Ta
%   and Tb with associated nodal coordinates in Na and Nb and
%   computes the interface projection matrix M

bl=[1];                        % bl: list of triangles of Tb to treat
bil=[1];                       % bil: list of triangles Ta to start with
bd=zeros(size(Tb,1)+1,1);      % bd: flag for triangles in Tb treated
bd(end)=1;                     % guard, to treat boundaries
bd(1)=1;                       % mark first triangle in b list
M=sparse(size(Nb,2),size(Na,2));
while length(bl)>0
  bc=bl(1); bl=bl(2:end);      % bc: current triangle of Tb
  al=bil(1); bil=bil(2:end);   % triangle of Ta to start with
  ad=zeros(size(Ta,1)+1,1);    % ad: flag for triangles in Ta treated
  ad(end)=1;                   % guard, to treat boundaries 
  ad(al)=1;                    % mark first triangle in a list
  n=[0 0 0];                   % triangles intersecting with neighbors
  while length(al)>0
    ac=al(1); al=al(2:end);    % take next candidate
    [P,nc,Mc]=Intersect(Nb(:,Tb(bc,1:3)),Na(:,Ta(ac,1:3)));
    if ~isempty(P)             % intersection found
      M(Tb(bc,1:3),Ta(ac,1:3))=M(Tb(bc,1:3),Ta(ac,1:3))+Mc;
      t=Ta(ac,3+find(ad(Ta(ac,4:6))==0));
      al=[al t];               % add neighbors
      ad(t)=1;                 % mark them as treated
      n(find(nc>0))=ac;        % ac is starting candidate for neighbor
    end
  end
  tmp=find(bd(Tb(bc,4:6))==0); % find non-treated neighbors
  idx=find(n(tmp)>0);          % take those which intersect with Ta
  t=Tb(bc,3+tmp(idx));
  bl=[bl t];                   % and add them
  bil=[bil n(tmp(idx))];       % with starting candidates Ta
  bd(t)=1;                     % mark them as treated 
end
