% QUADRAT A special quadrature routine used in R_MOD.
%
function xw=quadrat(N,i)
global ab0 M Z 
global theta c
xw=gauss(N,ab0);
if M>0
  for n=1:N
    p(n)=abs(prod((1+xw(n,1)*Z(1:M,1)).^Z(1:M,2)));
  end
  xw(:,2)=xw(:,2)./p';
%
% If an extra factor s is to be incorporated as in Remark to
% Theorem 3.25, then the preceding statement must be replaced by
%
% xw(:,2)=sqrt(1+.5*theta*xw(:,1)).*xw(:,2)./p';
%
end

