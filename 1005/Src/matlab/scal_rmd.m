% SCAL_RMD  Reverse mode differentiation of Blas operation SCAL
%   fx = SCAL_RMD(a,fx) updates fx for the operation x:=a*x (i.e. call
%   scal(a,x)). On entry fx contains the derivative of f w.r.t. the new x, and
%   on exit it contains the contribution to the derivative w.r.t. the old x
%   arising from the operation (if the old x has contributed to f in other ways,
%   fx should later be further updated).
%
%   [fx,fa] = SCAL_RMD(a,fx,x,fa) also adds to fa the current contribution to 
%   df/da (used when a is not constant).

function [fx,fa] = scal_rmd(a,fx,x,fa)
  if nargin > 2
    fa = fa + fx(:)'*x(:);  % fa = fa + dot(fx,x)
  end
  fx = a*fx(:);             % call scal(a,fx)
end
