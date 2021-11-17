function xw=gauss(N,ab)
%GAUSS Gauss quadrature rule.
%   Given a weight function w encoded by the Nx2 array AB of
%   the first N recurrence coefficients for the associated 
%   orthogonal polynomials, the first column of AB containing
%   the N alpha-coefficients and the second column the N beta-
%   coefficients, the call XW=GAUSS(N,AB) generates the nodes 
%   and weights XW of the N-point Gauss quadrature rule for the
%   weight function w. The nodes, in increasing order, are 
%   stored in the first column, the N corresponding weights in
%   the second column, of the Nx2 array XW.

N0=size(ab,1); if N0<N, error('input array ab too short'), end
J=zeros(N);
for n=1:N, J(n,n)=ab(n,1); end
for n=2:N
  J(n,n-1)=sqrt(ab(n,2));
  J(n-1,n)=J(n,n-1);
end
[V,D]=eig(J);
[D,I]=sort(diag(D));
V=V(:,I);
xw=[D ab(1,2)*V(1,:)'.^2];
