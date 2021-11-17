% ASUM_RMD  Reverse mode differentiation of Blas operation ASUM
%
%   fx = ASUMD(x,fx,fa) adds to fx the contribution to df/dx from the operation
%   a := asum(x) (sum of absolute values of elements of x). On entry fx contains
%   the already accumulated contribution to df/dx and fa contains df/da. If any
%   elements in x are zero the operation is not differentiable, and in this case
%   the corresponding elements of fx are not changed

function fx = asum_rmd(x,fx,fa)
  fx = fx(:) + fa*sign(x(:));
end
