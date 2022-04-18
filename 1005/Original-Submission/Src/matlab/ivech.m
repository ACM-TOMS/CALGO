% IVECH  inverse of vech operation
function L = ivech(x)
  n = fix(sqrt(length(x)*2));
  L = zeros(n);
  k = 1;
  for i=1:n
    k1 = k + n - i;
    L(i:n,i) = x(k:k+n-i);
    k = k1 + 1;
  end
end
  