%TEST_PRODUCT_DERIV  Test matrix product derivative functions
%
%  TEST_PRODUCT_DERIV checks that aat_deriv, ata_deriv and atb_deriv, which 
%  differentiate matrix products, work correctly.

function test_product_deriv
  fprintf('TESTING ATB_DERIV, ATA_DERIV AND AAT_DERIV...');
  A1 = rand(3,4);
  A2 = rand(3,4);
  B1 = rand(4,5);
  B2 = rand(4,5);
  x = [0.74; 0.33];
  d = diff_test(@fun, x, A1, A2, B1, B2, 2);
  ascertain(all(d<1e-8));
  disp('OK');
end

function [f,g] = fun(x,A1,A2,B1,B2,nPar)
  A = x(1)*A1 + x(2)*A2; Ad = cat(3, A1, A2); Adt = permute(Ad,[2,1,3]);
  B = x(1)*B1 + x(2)*B2; Bd = cat(3, B1, B2); Bdt = permute(Bd,[2,1,3]);
  [f{1},h{1}] = atb_deriv(A', Adt, B, Bd);    ascertain(almostequal(f{1}, A*B));
  [f{2},h{2}] = aat_deriv(A, Ad);             ascertain(Lequal(f{2}, A*A'));
  [f{3},h{3}] = ata_deriv(A, Ad);             ascertain(Lequal(f{3}, A'*A));  
  f{1} = f{1}(:);
  g{1} = reshape(h{1},length(f{1}),nPar);
  for i = 2:3
    f{i} = vech(f{i});
    for ip=1:nPar
      g{i}(:,ip) = vech(h{i}(:,:,ip));
    end
  end
end

function e = Lequal(A,B)
  e = almostequal(tril(A), tril(B));
end
