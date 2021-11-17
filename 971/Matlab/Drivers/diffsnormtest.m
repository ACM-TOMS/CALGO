%
% This script tests diffsnorm on dense matrices.
%
addpath('../src')

pass = [];

for m = [100 200]
  for n = [100 200]
    for isreal = [true false]

      if(isreal)
        A = randn(m,n);
      end
      if(~isreal)
        A = randn(m,n) + 1i*randn(m,n);
      end

      [U,S,V] = svd(A,'econ');
      snorm = diffsnorm(A,U,S,V,true);
      pass = [pass snorm<.1d-10*S(1,1)];

    end
  end
end


if(all(pass))
  disp('All tests succeeded.');
end

if(~all(pass))
  error('A test failed.');
end
