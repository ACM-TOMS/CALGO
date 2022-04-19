%
% This script tests adaptivepca on dense matrices.
%
addpath('../src')

err = [];
errs = [];

for m = [10 20]
  for n = [10 20]
    for k = [3 9]
      for b = [2 10]
        for p = [1 2]
          for isreal = [true false]

              if(isreal)
                U = randn(m,k);
                U = qr(U,0);
                V = randn(n,k);
                V = qr(V,0);
              end

              if(~isreal)
                U = randn(m,k) + 1i*randn(m,k);
                U = qr(U,0);
                V = randn(n,k) + 1i*randn(n,k);
                V = qr(V,0);
              end


              S0 = zeros(k,k);
              S0(1,1) = 1;
              S0(2,2) = .1;
              S0(3,3) = .01;


              A = U*S0*V';


              Ac = A;


              [U,S1,V] = svd(Ac,'econ');
              [U,S2,V] = adaptivepca(A,1.0d-13,b,p);
             

              S3 = zeros(min(m,n));
              S3(1:size(S2,1),1:size(S2,1)) = S2;

              errs = [errs norm(diag(S1)-diag(S3))];
              err = [err diffsnorm(A,U,S2,V,1)];
          end
        end
      end
    end
  end
end

errs;
err;


if(all(err<.1d-10))
  disp('All tests succeeded.');
end

if(~all(err<.1d-10))
  error('A test failed.');
end
