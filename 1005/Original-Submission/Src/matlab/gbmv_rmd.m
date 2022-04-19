% GBMV_RMD  Reverse mode differentiation of Blas-2 operation gbmv
%
% fx = gbmv_rmd(Ab, x, kl, alpha, beta, trans, fy, fx) adds to fx the
% contribution to df/dx from the operation y = alpha*op(A)*x + beta*y (i.e. the
% Blas call gbmv(A,x,y,kl,m,alpha,beta,trans)). A is a band matrix stored in Ab
% in Blas fashion (see below). On entry fy is df/dy (differentiated w.r.t. the
% new y) and  fx contains the contributions to df/dx already accumlated.
%
% [fx,fy] = gbmv_rmd(Ab, x, kl, alpha, beta, trans, fy, fx) sets in addition fy
% to contain the contribution to the derivative w.r.t. the old y arising from
% the operation (if the old y has contributed to f in other ways, fy should
% later be further updated).
%
% [fx,fy,fA] = gbmv_rmd(Ab, x, kl, alpha, beta, trans, fy, fx, fA) also adds to
% fA the contribution to df/dA from the operation.
%
% The rows of Ab store the diagonals of A from top to bottom; kl is the number
% lower diagonals and the number of upper diagonals is ku = size(A,1) - kl - 1).
% The upper diagonals are stored right-aligned and the lower ones are
% left-aligned (so the upper left ku by ku triangle of Ab is unused, and,
% depending on dimensions, possibly a lower right triangle also).
%
% If alpha and/or beta are 0 or 1 these facts are taken advantage of.
% [fx,dummy,fA] = gbmv_rmd(Ab, x, kl, alpha, 0, trans, 0, fx, fA) can be used
% for beta = 0.

function [fx, fy, fA] = gbmv_rmd(Ab, x, kl, alpha, beta, trans, fx, fy, fA)
  x = x(:);
  fx = fx(:);
  fy = fy(:);
  if trans=='T', transt = 'N'; else transt = 'T'; end
  if alpha == 1
    fx = gbmv(Ab, fy, fx, kl, 1, 1, transt);  % fx = fx + op(A)'*fy
    if exist('fA', 'var')
      fA = gbr(fA, fy, x, kl, 1);             % fA = fA + fy*x'
    end
  elseif alpha ~= 0
    fx = gbmv(Ab, fy, fx, kl, alpha, 1, transt);  % fx = fx + op(A)'*fy
    if exist('fA', 'var')
      fA = gbr(fA, fy, x, kl, alpha);             % fA = fA + fy*x'
    end
  end
  if beta == 0 && nargout > 1
    fy = zeros(size(fy));             % fy := [0,0,...,0]'
  elseif beta ~= 1 && nargout > 1
    fy = beta*fy;                     % call scal(beta,fy)
  end
end
  
function y = gbmv(Ab, x, y, kl, alpha, beta, trans)
  % Sets y := alpha*op(A)*x + beta*y where op(A) is A' if trans is 'T',
  % otherwise op(A) = A. A is banded, stored in Ab as described above.
  [kd,n] = size(Ab);
  ku = kd - 1 - kl;
  if beta~=1, y = y*beta; end
  if alpha==0, return, end
  if alpha~=1, x = x*alpha; end
  if trans=='T'
    m = length(x);
    assert(length(y)==n)
    for j=1:m
      i1 = max(1,j-ku);
      i2 = min(n,j+kl);
      k1 = max(1,ku+2-j);
      k2 = k1 + (i2 - i1);
      y(i1:i2) = y(i1:i2) + x(j)*Ab(j,k1:k2)';
    end
  else
    m = length(y);
    assert(length(x)==n)
    for j=1:n
      i1 = max(1,j-ku);
      i2 = min(m,j+kl);
      k1 = max(1,ku+2-j);
      k2 = k1 + (i2 - i1);
      y(i1:i2) = y(i1:i2) + x(j)*Ab(k1:k2,j);
    end
  end
end

function  Ab = gbr(Ab, x, y, kl, alpha)
  % A special version of Blas _ger subroutine, that updates a banded matrix with
  % the corresponding part of x*y'.
  [kd,n] = size(Ab);
  ku = kd - 1 - kl;
  m = length(x);
  assert(length(y)==n);
  if alpha==0, return, end
  if alpha~=1, x = x*alpha; end
  for j=1:n
    i1 = max(1,j-ku);
    i2 = min(m,j+kl);
    k1 = max(1,ku+2-j);
    k2 = k1 + (i2 - i1);
    Ab(k1:k2,j) = Ab(k1:k2,j) + x(i1:i2)*y(j);
  end
end
