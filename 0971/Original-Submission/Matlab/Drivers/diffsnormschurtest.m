%
% This script tests diffsnormschur on dense matrices.
%
addpath('../src')

pass = [];

for n = [100 200]
  for isreal = [true false]

    if(isreal)
      A = randn(n,n);
    end
    if(~isreal)
      A = randn(n,n) + 1i*randn(n,n);
    end

    [U,S,V] = svd(A,'econ');
    T = S*V'*U;
    snorm = diffsnormschur(A,U,T);
    pass = [pass snorm<.1d-10*S(1,1)];

  end
end


if(all(pass))
  disp('All tests succeeded.');
end

if(~all(pass))
  error('A test failed.');
end
