% GEMV_RMD  Reverse mode differentiation of Blas-2 operation GEMV
%
% xa = gemv_rmd(alpha, A, x, beta, xa, ya) adds to the adjoint vector xa
% the contribution to df/dx from the operation y = alpha*A*x + beta*y (i.e.
% call gemv(alpha,A,x,beta,y)). On entry ya is df/dy (differentiated w.r.t.
% the new y) and xa contains the contributions to df/dx already accumlated.
%
% [xa,ya] = gemv_rmd(alpha, A, x, beta, xa, ya) sets in addition ya to
% contain the contribution to the adjoint of y arising from the operation
% (if the old y has contributed to f in other ways, ya should later be
% further updated).
%
% [xa,ya,Aa] = gemv_rmd(alpha, A, x, beta, xa, ya, Aa) also adds to Aa the
% contribution to df/dA from the operation.
%
% If the final output result f has more than one element then the adjoint
% on each one is updated (all of xa, ya and Aa have an additional
% (rightmost) dimension which is looped over).
%
% If alpha and/or beta are 0 or 1 these facts are taken advantage of.
% [xa,dummy,Aa] = gemv_rmd(alpha, A, x, 0, xa, 0, Aa) can be used for beta
% = 0.

function [xa, ya, Aa] = gemv_rmd(a, A, x, b, xa, ya, Aa)
  [~,~,nf] = size(Aa);
  if a == 1                                % CORRESPONDING FORTRAN 95 CALLS:
    xa = xa + A'*ya;                       % call gemv(1,A',ya,1,xa)
    if nargin > 6                          % or gemm(...)
      for j = 1:nf
        Aa(:,:,j) = Aa(:,:,j) + ya(:,j)*x';% call ger(1,ya(.,j),x,Aa(.,.,j)
      end
    end
  elseif a ~= 0
    xa = xa + a*(A'*ya);                   % call gemv(a,A',ya,1,xa)
    if nargin > 6
      for j = 1:nf
        Aa(:,:,j) = Aa(:,:,j) + a*ya(:,j)*x'; % call ger(a,ya(.,j),x,Aa(.,.,j)
      end
    end
  end
  if b == 0 && nargout > 1
    ya = zeros(size(ya));                  % ya := [0,0,...,0]'
  elseif b ~= 1 && nargout > 1
    ya = b*ya;                             % call scal(b,ya)
  end
end
