%
% This script tests adaptivepca on sparse matrices.
%
addpath('../src')

err = [];
errs = [];
bests = [];

for m = [100 200]
  for n = [100 200]
    for k = [30 90]
      for b = [5 10]
        for P = [1 2]
          for isreal = [true false]

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


              Ac = A;


              [U,S1,V] = svd(full(Ac),'econ');
              [U,S2,V] = adaptivepca(A,1.0d-10,b,P);


              S3 = zeros(min(m,n));
              S3(1:size(S2,1),1:size(S2,1)) = S2;
              errs = [errs norm(diag(S1)-diag(S3))];


              err = [err diffsnorm(A,U,S2,V,1)];


              bests = [bests S1(k+1,k+1)];
          end
        end
      end
    end
  end
end



if(all(err./bests<10 | err<1d-10))
  disp('All tests succeeded.');
end

if(~all(err./bests<10 | err<1d-10))
  error('A test failed.');
end
