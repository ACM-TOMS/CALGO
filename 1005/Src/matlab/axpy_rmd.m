% AXPY_RMD  Reverse mode differentiation of Blas operation AXPY
%   fx = AXPYD(a,x,fx,fy) adds to fx the contribution to df/dx from the Blas
%   operation y = a*x + y (i.e. call axpy(a,x,y)). On entry fy is df/dy and fx
%   contains the contributions to df/dx already accumlated.
%
%   [fx,fa] = AXPYD(a,x,fx,fy,fa) also updates df/da (used when a is not
%   constant).

function [fx,fa] = axpy_rmd(a,x,fy,fx,fa)
  fy = fy(:);
  fx = a*fy + fx(:);     % call axpy(a,fy,fx)
  if nargin > 4
    fa = fa + fy'*x(:);  % fa = fa + dot(fy,x)
  end
end
