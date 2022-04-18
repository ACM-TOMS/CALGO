function InterfaceMatrix3d(Na,Ta,Nb,Tb);
% INTERFACEMATRIX3D projection matrix for nonmatching thetrahedra grids 
%   M=InterfaceMatrix3d(Na,Ta,Nb,Tb); takes two tetrahedra meshes Ta
%   and Tb with associated nodal coordinates in Na and Nb and
%   computes the projection matrix M. As a precondition, the first
%   thetrahedra in Ta and Tb need to intersect.

bl=[1];                        % bl: list of tetrahedra of Tb to treat
bil=[1];                       % bil: list of tetrahedra Ta to start with
bd=zeros(size(Tb,1)+1,1);      % bd: flag for tetrahedra in Tb treated 
bd(end)=1;                     % guard, to treat boundaries
bd(1)=1;                       % mark first tetrahedra in b list.
while length(bl)>0
  bc=bl(1); bl=bl(2:end);      % bc: current tetrahedra of Tb 
  al=bil(1); bil=bil(2:end);   % tetrahedra of Ta to start with
  ad=zeros(size(Ta,1)+1,1);    % same as for bd
  ad(end)=1;
  ad(al)=1; 
  n=[0 0 0 0];                 % tetrahedra intersecting with neighbors
  while length(al)>0
    ac=al(1); al=al(2:end);    % take next candidate
    [P,nc,H,v]=Intersection3d(Nb(:,Tb(bc,1:4)),Na(:,Ta(ac,1:4)));
    if ~isempty(P)             % intersection found
      t=Ta(ac,4+find(ad(Ta(ac,5:8))==0));
      al=[al t];               % add neighbors
      ad(t)=1;
      n(find(nc>0))=ac;        % ac is starting candidate for neighbor  
    end
  end
  tmp=find(bd(Tb(bc,5:8))==0); % find non-treated neighbors
  idx=find(n(tmp)>0);          % only if intersecting with Ta
  t=Tb(bc,4+tmp(idx));
  bl=[bl t];                   % and add them
  bil=[bil n(tmp(idx))];       % with starting candidates Ta
  bd(t)=1;
end
