%
% This script tests diffsnormschur on sparse matrices.
%
addpath('../src')

pass = [];

for n = [100 200]
  for isreal = [true false]

    if(isreal)
      A = 2*spdiags((1:n)',0,n,n);
    end
    if(~isreal)
      A = 2*spdiags((1:n)',0,n,n)*(1+1i);
    end

    A = A - spdiags((0:(n-1))',1,n,n);
    A = A - spdiags((1:n)',-1,n,n);

    [U,S,V] = svd(full(A),'econ');
    A = A/S(1,1);

    [U,S,V] = svd(full(A),'econ');
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
