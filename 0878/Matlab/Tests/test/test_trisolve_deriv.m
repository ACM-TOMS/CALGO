%TEST_TRISOLVE_DERIV  Test functions forward_sub_deriv and back_sub_deriv
%
%  TEST_TRISOLVE_DERIV does some preliminary testing of the trianglar solve
%  derivative functions forward_sub_deriv and back_sub_deriv. Further tests are
%  done by e.g. test_chol_deriv.

function test_trisolve_deriv
  fprintf('PRELIMINARY TEST OF FORWARD_SUB_DERIV AND BACK_SUB_DERIV...')
  L = 1;
  Ld = ones(1,1,2);
  Z = 1;
  Zd = ones(1,1,2);
  X = L\Z;
  Y = L'\Z;
  Xd = forward_sub_deriv(L,Ld,X,Zd);
  Yd = back_sub_deriv(L,Ld,Y,Zd);
  ascertain(norm(Xd(:,:,1))<1e-14);
  ascertain(norm(Xd(:,:,2))<1e-14);
  ascertain(norm(Yd(:,:,1))<1e-14);
  ascertain(norm(Yd(:,:,2))<1e-14);
  %
  L = chol(hilb(2))';
  Ld = zeros(2,2,3);
  for i=1:3, Ld(:,:,i) = tril(rand(2,2)); end
  Z = rand(2,2);
  Zd = rand(2,2,3);
  X = L\Z;
  Y = L'\Z;
  Xd = forward_sub_deriv(L,Ld,X,Zd);
  Yd = back_sub_deriv(L,Ld,Y,Zd);
  for i=1:3
    ascertain(almostequal(L*Xd(:,:,i) + Ld(:,:,i)*X, Zd(:,:,i)));
    ascertain(almostequal(L'*Yd(:,:,i) + Ld(:,:,i)'*Y, Zd(:,:,i)));
  end
  disp('OK')
end
