%
% This script tests pcafast on sparse matrices.
%
addpath('../src')

err = [];
errs = [];
bests = [];

for m = [100 200]
  for n = [100 200]
    for k = [30 90]
      for its = [2 1000]
        for raw = [true false]
          for isreal = [true false]


            l = k+1;


            if(isreal)
              A = 2*spdiags((1:min(m,n))',0,m,n);
            end
            if(~isreal)
              A = 2*spdiags((1:min(m,n))',0,m,n)*(1+1i);
            end


            A = A - spdiags((0:(min(m,n)-1))',1,m,n);
            A = A - spdiags((1:min(m,n))',-1,m,n);
            A = A/sqrt(normest(A*A'));
            A = A*A'*A;
            A = A*A'*A;
            A(randperm(m),:) = A;
            A(:,randperm(n)) = A;


            if(raw)
              Ac = A;
            end
            if(~raw)
              c = sum(A)/m;
              Ac = A - ones(m,1)*c;
            end


            [U,S1,V] = svd(full(Ac),'econ');
            [U,S2,V] = pcafast(A,k,raw,its,l);


            S3 = zeros(min(m,n));
            S3(1:k,1:k) = S2;
            errs = [errs norm(diag(S1)-diag(S3))];


            err = [err diffsnorm(A,U,S2,V,raw)];


            bests = [bests S1(k+1,k+1)];


          end
        end
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
