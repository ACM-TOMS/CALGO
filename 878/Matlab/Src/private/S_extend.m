% S_EXTEND  Extend S0...Sp from Yule-Walker equations to S0...S(n-1)
%
%   Scol = S_EXTEND(A, G, S, n) when A = [A0...Ap], G = {G0,...,Gq} and S =
%   {S0,...,Sp} calculates S(p+1)...S(n-1) and returns them together with
%   S0,..., Sp in a block column vector Scol = [S0; S1;...; S(n-1)]. Thus Scol
%   will be n·r × r.
%
%   [Scol, Scold] = S_EXTEND(A, G, S, n, Gd, Sd) returns also the n·r × r × nPar
%   derivatives of Scol.

function [Scol,Scold] = S_extend(A, G, S, n, Gd, Sd)
  DIFF = nargout > 1;
  r = size(G{1}, 1);
  p = size(A,2)/r;
  q = length(G)-1;
  A = cell2mat(fliplr(makecell(A))); % change A from [A1..Ap] to [Ap..A1]
  if isempty(A), A = zeros(r,0); end
  Scol = cat(1, S{:}, zeros(r*(n-p-1),r));
  if DIFF
    nPar = size(Sd{1},3);
    Scold = cat(1, Sd{:}, zeros(r*(n-p-1),r,nPar));
  end
  pp = p+1;
  K = r+1:r*pp;
  I = r*pp+1:r*(pp+1);
  for j = p+1:n-1
    if j <= q
      Scol(I,:) = G{j+1}; 
      if DIFF, Scold(I,:,:) = Gd{j+1}; end
    end
    Scol(I,:) = Scol(I,:) + A*Scol(K,:);
    if DIFF, Scold(I,:,:) = Scold(I,:,:) + AGdiff(A, Scol(K,:), Scold(K,:,:)); end
    K = K+r;
    I = I+r;
  end
end

function Fd = AGdiff(A, G, Gd)
  nPar = size(Gd, 3);
  r = size(A,1);
  p = size(A,2)/r;
  Fd = reshape(A*reshape(Gd,r*p,r*nPar),r,r,nPar);
  k = 1;
  for i = p-1:-1:0
    for c=r*i+1:r*i+r
      for l=1:r
        Fd(l,:,k) = Fd(l,:,k) + G(c,:);
        k = k+1;
      end
    end
  end
end
