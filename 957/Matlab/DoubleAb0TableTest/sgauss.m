% SGAUSS Symbolic counterpart of gauss.m
%   The array ab must have dimension Nx2. 
%   Otherwise, if the dimension is N1x2
%   with N1>N, then the routine should be
%   called with ab(1:N,:) instead of ab.
%
function xw=sgauss(dig,N,ab)
N0=size(ab,1); if N0<N, error('input array ab too short'), end
sab=vpa(ab,dig);
if N==1
  xw=[sab(1,1) sab(1,2)];
  return
elseif N==2
  J=[sab(1,1) sqrt(sab(2,2)); sqrt(sab(2,2)) sab(2,1)];
else
  J=diag(sqrt(sab(2:end,2)),1);
  J=diag(sab(:,1))+J+J';
end
[V,D]=eig(J);
D=diag(D);
[v,I]=sort(double(D));
D=D(I);
V=V(:,I);
xw=[D vpa(sab(1,2))*V(1,:)'.^2];
