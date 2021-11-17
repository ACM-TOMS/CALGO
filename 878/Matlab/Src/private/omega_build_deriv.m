% OMEGA_BUILD_DERIV  Build derivatives of matrix Omega
% 
%  [Sud,Olowd] = OMEGA_BUILD_DERIV(Sd,Gd,Wd,n) returns an h·r×h·r×nPar array Sud
%  and an (n-h)r×(q+1)r×nPar array Olowd with the derivatives of Su and Olow
%  with respect to all parameters. Sd, Gd and Wd are the derivatives of the
%  Si, Gi and Wi matrices.

function [Sud, Olowd] = omega_build_deriv(Sd,Gd,Wd,p,n)
  q = length(Gd) - 1;
  h = max(p,q);
  r = size(Gd{1},1);
  Sd = cell2mat(fliplr(Sd(1:p)));
  Gd = cell2mat(fliplr(Gd));
  Wd = cell2mat(fliplr(Wd));
  nPar = size(Gd,3);  
  Sud = zeros(r*h, r*h, nPar);
  Olowd = zeros(r*(n-h), r*(q+1), nPar);
  K = 1:r;
  for t = 1:p
    Sud(K, 1:t*r, :) = Sd(:, end-t*r+1:end, :);
    K = K + r;
  end
  J = (q-p)*r+1 : q*r;
  for t = p+1:h
    Sud(K, 1:p*r, :) = Gd(:, J, :);
    Sud(K, p*r+1:t*r, :) = Wd(:, end-(t-p)*r+1:end, :);
    J = J - r;
    K = K + r;
  end
  K = 1:r;
  je = min(p,q)*r;
  for t = 1:n-h
    Olowd(K, 1:je, :) = Gd(:, 1:je, :);
    Olowd(K, je+1:end, :) = Wd(:, je+1:end, :);
    if je > 0, je = je-r; end
    K = K + r;
  end
end
