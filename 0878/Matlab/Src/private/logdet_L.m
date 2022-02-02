%LOGDET_L  Logarithm of determinant of a lower triangular matrix
%
%  ld = LOGDET_L(L) returns log(det(L)) and [ld, ldd] = LOGDET_L(L, Ld) finds in
%  addition the derivative of ld with respect to several parameters. If L is
%  n×n, Ld should be n×n×nPar and contain the derivative of L (nPar is the
%  number of parameters)

function [ld, ldd] = logdet_L(L, Ld);
  ld = sum(log(diag(L)));
  if nargin>1
    nPar = size(Ld,3);
    ldd = zeros(1,nPar);
    for j=1:nPar
      ldd(j) = sum(diag(Ld(:,:,j))./diag(L));
    end
  end
end
