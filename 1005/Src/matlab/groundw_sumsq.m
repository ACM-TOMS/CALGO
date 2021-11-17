function [theta, ka] = groundw_sumsq(k, f0, F, A0, A, C, d)
  % k is an m-vector of parameters to be determined to minimize theta
  %
  % The remaining parameters are data:
  %   f0 is an n-vector
  %   F is n by m
  %   A0 is n by n
  %   A{i} is n by n for i=1...m
  %   C is s by n
  %   d is s by 1
  
  % COMPUTE LIKELIHOOD
  m = length(k);
  B = A0;
  for j = 1:m, B = B + k(j)*A{j}; end   % axpy(k(j), vec(A{j}), vec(B))
  f = F*k + f0;                         % gemv
  [L,U,p] = lu(B, 'vector');
  u = U\(L\f(p));                       % solve linear system B*u = f
  r = C*u - d;                          % gemv
  theta = r'*r;                         % dot
  
  if nargout > 1
    % COMPUTE ADJOINT OF k
    thetaa = 1;
    ra = 2*thetaa*r;                      % dot_rmd
    ua = C'*ra;                           % gemv_rmd
    z(p) = L'\(U'\ua); z = z(:);          % linsys_rmd (solve B'*z = ua)
    fa = z;                               % linsys_rmd
    Ba = -z*u';                           % linsys_rmd (ger)
    ka = F'*fa;                           % gemv_rmd
    for j=m:-1:1
      ka(j) = ka(j) + dot(Ba(:),A{j}(:)); % axpy_rmd
    end
  end
end
