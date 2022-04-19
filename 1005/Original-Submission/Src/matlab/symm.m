% SYMM  Symmetric matrix with given L part
%
%   A = SYMM(L) for matrix L returns the symmetric matrix with lower part L
%   A = SYMM(x) for vector x returns symm(ivech(x))
%
% Note: Had to rename sym to symm because of a but in Matlab R2018a

function A = symm(L)
  if isvector(L), L = ivech(L); end
  A = tril(L) + tril(L)' - dg(L);
end
