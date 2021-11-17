% ATBA_C  Calculate A'·B·A + C and optionally its derivatives
%
%   D = ATBA_C(A,B,C) calculates D = A'·B·A + C. B and C are assumed to be
%   symmetric, only their lower triangle is used (and only the lower triangle of
%   D is guaranteed to be correct).
%
%   [D,Dd] = ATBA_C(A,B,C,Ad,Bd,Cd) calculates also the derivative w.r.t.
%   parameter theta(i) in Dd(:,:,i), i = 1,...,nPar. Ad(:,:,i), Bd(:,:,i) and
%   Cd(:,:,i) should be the derivatives of A, B and C w.r.t. theta(i).

function [D, Dd] = atba_c(A, B, C, Ad, Bd, Cd)
  DIFF = nargin > 3;
  if ~DIFF
    L = tril(B);
    for i=1:size(B,1), L(i,i) = L(i,i)/2; end
    ALA = A'*L*A;
    D = ALA + ALA' + C;
  else
    [m,n,nPar] = size(Bd);
    L = tril(B);
    Ld = Bd;
    for i=1:m
      L(i,i,:) = L(i,i,:)/2;
      Ld(i,i,:) = Ld(i,i,:)/2; 
      Ld(1:i-1,i,:) = 0; 
    end
    [LA, LAd] = atb_deriv(L,Ld,A,Ad);
    [ALA,ALAd] = atb_deriv(A,Ad,LA,LAd);
    D = ALA + ALA' + C;
    Dd = Cd + ALAd + permute(ALAd,[2,1,3]);
  end
end
