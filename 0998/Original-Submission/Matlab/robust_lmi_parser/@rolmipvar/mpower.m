function Z = mpower(X, n)
%MPOWER (overloaded)
%
% Author: Alexandre Felipe
% 2014, Dec, 9

if( floor(n) ~= n )
  error('Unable to raise rolmipvar to non-integer power');
elseif (n < 0)
  error('Exponent must be a positive integer.');
else
  Z = 1;
  for i = 2:n
    Z = mtimes(X, Z);
  end
end
end

