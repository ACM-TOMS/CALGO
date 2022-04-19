function [x,y] = TwoProduct(a,b)
%TWOPRODUCT   Error free transformation of a*b into x+y with x=fl(a*b)
%
%   [x,y] = TwoProduct(a,b)
%
%On return, x+y=a*b and x=fl(a*b) provided no over- or underflow occurs .
%Input a,b may be vectors or matrices as well, in single or double precision.
%
%Follows G.W. Veltkamp, see T.J. Dekker: A floating-point technique for 
%  extending the available precision, Numerische Mathematik 18:224-242, 1971.
%Requires 17 flops for scalar input.
%
% written  03/03/07     S.M. Rump
%

    x = a.*b;
    if any(~isfinite(x(:))) 
      error('overflow occurred in TwoProduct')
    end
    [ah,al] = Split(a);
    [bh,bl] = Split(b);
    y = al.*bl - ( ( ( x - ah.*bh ) - al.*bh ) - ah.*bl );
end
