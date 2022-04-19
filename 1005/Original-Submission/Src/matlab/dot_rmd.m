% DOT_RMD  Reverse mode differentiation of Blas operation DOT
%
%   [fx,fy] = DOT_RMD(x,y,fx,fy,fa) adds to fx and fy the contribution to df/dx
%   and df/dy from the Blas operation a = x'*y (i.e. a = dot(x,y)). On entry fa
%   is df/da, and fx and fy contain the contributions to df/dx already
%   accumulated. Any input parameter row vectors are changed to column vectors.
%
%   fx = dot_rmd(x,fx,fa) updates fx with the contribution from the Blas
%   operation a = x'*x.

function [fx,fy] = dot_rmd(x,varargin)
  if nargin == 2
    [fx,fa] = deal(varargin{:});
    fx = fx(:) + (2*fa)*x(:);          % call axpy(2*fa, x, fx)
  else
    [y,fx,fy,fa] = deal(varargin{:});
    fx = fx(:) + fa*y(:);              % call axpy(fa, y, fx)
    fy = fy(:) + fa*x(:);              % call axpy(fa, x, fy)    
  end
end
