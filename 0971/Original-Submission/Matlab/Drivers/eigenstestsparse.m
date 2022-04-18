%
% This script tests eigen on sparse self-adjoint matrices.
%
addpath('../src')

err = [];
errs = [];
bests = [];

for n = [100 200]
  for k = [30 90]
    for its = [2 1000]
      for isreal = [true false]


        l = k+1;


        A = 2*spdiags((1:(n/2))',0,n/2,n/2);


        if(isreal)
          A = A - spdiags((0:((n/2)-1))',1,n/2,n/2);
          A = A - spdiags((1:(n/2))',-1,n/2,n/2);
        end
        if(~isreal)
          A = A - 1i*spdiags((0:((n/2)-1))',1,n/2,n/2);
          A = A + 1i*spdiags((1:(n/2))',-1,n/2,n/2);
        end

        A = A/normest(A);
        A = A*A;
        A = A*A;
        A = A*A;
        A = A*A;
        A = [zeros(n/2,n/2) A; A zeros(n/2,n/2)];
        P = randperm(n);
        A(P,:) = A;
        A(:,P) = A;


        [V,D1] = eig(full(A));
        [V,D2] = eigen(A,k,its,l);


        D3 = zeros(n);
        D3(1:k,1:k) = D2;

        [E,P] = sort(abs(diag(D1)),'descend');
        D1 = diag(D1);
        D1 = diag(real(D1(P)));

        errs = [errs norm(abs(diag(D1))-abs(diag(D3)))];


        err = [err diffsnormschur(A,V,D2)];


        bests = [bests abs(D1(k+1,k+1))];


      end
    end
  end
end


if(all(err./bests<10 | err<.1d-10))
  disp('All tests succeeded.');
end

if(~all(err./bests<10 | err<.1d-10))
  error('A test failed.');
end
