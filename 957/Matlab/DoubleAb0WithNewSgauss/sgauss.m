function xw=sgauss(dig,N,ab)
%SGAUSS Symbolic dig-digit precision counterpart of gauss.m.

digits(dig)
N0=size(ab,1); if N0<N, error('input array ab too short'), end
sab=vpa(ab);
J=zeros(N); J=vpa(J);
for n=1:N, J(n,n)=sab(n,1); end
for n=2:N
  J(n,n-1)=sqrt(sab(n,2));
  J(n-1,n)=J(n,n-1);
end
[V,D]=eig(J);
D=diag(D);
[v,I]=sort(double(D));
D=D(I);
V=V(:,I);
xw=[D vpa(sab(1,2))*V(1,:)'.^2];
