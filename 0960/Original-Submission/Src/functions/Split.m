function [x,y] = Split(a)
%SPLIT        Error-free split a=x+y into two parts.
%
%   [x,y] = Split(a)
%
%On return, x+y=a and both x and y need at most k bits in the mantissa.
%In double precision k=26, in single precision k=12.
%Input may be a vector or matrix as well.
%
%Follows T.J. Dekker: A floating-point technique for extending the available
%  precision, Numerische Mathematik 18:224-242, 1971.
%Requires 4 flops for scalar input.
%
%
% written  03/03/07     S.M. Rump
%

  if isa(a,'double'), prec='double'; else prec='single'; end    
  factor = 2^ceil(log2(2/eps(prec))/2)+1;

  c = factor*a;            % factor('double')=2^27+1, factor('single')=2^12+1
  if any(~isfinite(c(:)))
    error('overflow in Split')
  end
  x = c - ( c - a );
  y = a - x;
