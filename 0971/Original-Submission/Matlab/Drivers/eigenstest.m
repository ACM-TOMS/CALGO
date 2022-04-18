%
% This script tests eigen on dense self-adjoint matrices.
%
addpath('../src')

err = [];
errs = [];

for n = [10 20]
  for k = [3 9]
    for its = [0 2 1000]
      for isreal = [true false]


        l = k+1;


        if(isreal)
          V = randn(n,k);
        end
        if(~isreal)
          V = randn(n,k) + 1i*randn(n,k);
        end

        V = qr(V,0);


        D0 = zeros(k,k);
        D0(1,1) = 1;
        D0(2,2) = -.1;
        D0(3,3) = .01;


        A = V*D0*V';
        A = (A+A')/2;


        [V,D1] = eig(A);
        [V,D2] = eigen(A,k,its,l);


        D3 = zeros(n);
        D3(1:k,1:k) = D2;

        [E,P] = sort(abs(diag(D1)),'descend');
        D1 = diag(D1);
        D1 = diag(real(D1(P)));

        errs = [errs norm(diag(D1)-diag(D3))];


        err = [err diffsnormschur(A,V,D2)];


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
