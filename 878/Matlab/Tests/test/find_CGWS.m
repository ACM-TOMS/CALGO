% FIND_CGWS
%
%   [C,G,W,S] = find_CGWS(A, B, Sig) finds the matrices Cj, Gj, Wj and Sj in
%   cell arrays {C1...Cq}, {G1...Gq}, {W1...Wq} and {S1...Sp}.
%
%   [C,G,W,S,Cd,Gd,Wd,Sd] = find_CGWS(A, B, Sig) returns also derivatives
%
%   The routine is used by the testing suite. It calls find_CGW (and
%   find_CGW_deriv) and the vyw-functions. See further comments in these.

function [C,G,W,S,Cd,Gd,Wd,Sd] = find_CGWS(A, B, Sig);
  if nargout <=4
    [C,G,W] = find_CGW(A, B, Sig);
    LUvyw = vyw_factorize(A);
    S = vyw_solve(A, LUvyw, G);
  else
    [CCd, GGd, WWd] = find_CGW_deriv(A, B, Sig);
    [C, Cd] = der2array(CCd);
    [G, Gd] = der2array(GGd);
    [W, Wd] = der2array(WWd);
    LUvyw = vyw_factorize(A);
    S = vyw_solve(A, LUvyw, G);
    RHS = vyw_deriv_rhs(A, GGd, S);
    Sd = vyw_solve(A, LUvyw, RHS);
  end
end
