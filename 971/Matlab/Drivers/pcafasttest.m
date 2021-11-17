%
% This script tests pcafast on dense matrices.
%
addpath('../src')

err = [];
errs = [];

for m = [10 20]
  for n = [10 20]
    for k = [3 9]
      for its = [0 2 1000]
        for raw = [true false]
          for isreal = [true false]


            l = k+1;


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


            if(raw)
              Ac = A;
            end
            if(~raw)
              c = sum(A)/m;
              Ac = A - ones(m,1)*c;
            end


            [U,S1,V] = svd(Ac,'econ');
            [U,S2,V] = pcafast(A,k,raw,its,l);


            S3 = zeros(min(m,n));
            S3(1:k,1:k) = S2;
            errs = [errs norm(diag(S1)-diag(S3))];


            err = [err diffsnorm(A,U,S2,V,raw)];


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
