%TEST_FIND_CGW  Test find_CGW which finds lagged cross-covariance matrices
%
%  TEST_FIND_CGW prints out "OK" if find_CGW works OK for three test cases with
%  three different values of (p,q,r) (find_CGW finds the (lagged)
%  cross-covariance matrices of (xt,epst), (yt,xt) and (yt,yt)). It may also be
%  useful for debugging. TEST_FIND_CGW also does a partial test the function
%  find_CGW_deriv, to be exact it tests that the function value returned by it
%  is the same as find_CGW returns. Further tests are carried out by
%  test_vyw_deriv.
%
function test_find_CGW
  fprintf('TESTING FIND_CGW... ')
  A1 = [1 1; 1 1];
  A2 = [2 2; 1 1];
  B1 = [1 1; 2 2];
  B2 = [2 2; 2 2];
  B3 = [3 3; 2 2];
  Sig = [3 1; 1 3];
  r = size(A1,1);
  I = eye(r);
  %
  % check p=2, q=3:
  [C,G,W] = find_CGW([A1,A2],[B1,B2,B3],Sig);
  [CCd,GGd,WWd] = find_CGW_deriv([A1,A2],[B1,B2,B3],Sig);
  ascertain(almostequal({CCd.mat}, C) && almostequal({GGd.mat}, G) && almostequal({WWd.mat},W));
  %
  Wb = W;
  C0 = C{1}; C = C(2:end);
  G0 = G{1}; G = G(2:end);  
  W0 = W{1}; W = W(2:end);  
  ascertain(length(G) == 3);
  ascertain(length(W) == 3);
  ascertain(almostequal(G0,     Sig + B1*C{1}' + B2*C{2}' + B3*C{3}'));
  ascertain(almostequal(G{1},B1*Sig + B2*C{1}' + B3*C{2}'           ));
  ascertain(almostequal(G{2},B2*Sig + B3*C{1}'                      ));
  ascertain(almostequal(G{3},B3*Sig                                 ));
  B0 = I;
  ascertain(almostequal(W0  , B0*Sig*B0' + B1*Sig*B1' + B2*Sig*B2' + B3*Sig*B3'));
  ascertain(almostequal(W{1}, B1*Sig*B0' + B2*Sig*B1' + B3*Sig*B2'             ));
  ascertain(almostequal(W{2}, B2*Sig*B0' + B3*Sig*B1'                          ));
  ascertain(almostequal(W{3}, B3*Sig*B0'                                       ));
  %
  % check p = 0, q=3
  [C,G,W] = find_CGW([],[B1,B2,B3],Sig);
  [CCd,GGd,WWd] = find_CGW_deriv([],[B1,B2,B3],Sig);
  ascertain(almostequal({CCd.mat}, C) && almostequal({GGd.mat}, G) && almostequal({WWd.mat},W));
  %
  C0 = C{1}; C = C(2:end);
  G0 = G{1}; G = G(2:end);
  ascertain(almostequal(C0,Sig));
  ascertain(length(G) == 3);
  ascertain(almostequal(W,Wb));
  ascertain(almostequal(G0  ,   Sig + B1*C{1}' + B2*C{2}' + B3*C{3}'));
  ascertain(almostequal(G{1},B1*Sig + B2*C{1}' + B3*C{2}'           ));
  ascertain(almostequal(G{2},B2*Sig + B3*C{1}'                      ));
  ascertain(almostequal(G{3},B3*Sig                                 ));
  %
  % check p=2, q=0:
  [C,G,W] = find_CGW([A1,A2],[],Sig);
  [CCd,GGd,WWd] = find_CGW_deriv([A1,A2],[],Sig);
  ascertain(almostequal({CCd.mat}, C) && almostequal({GGd.mat}, G) && almostequal({WWd.mat},W));
  %
  C0 = C{1}; C = C(2:end);
  G0 = G{1}; G = G(2:end);  
  W0 = W{1}; W = W(2:end);  
  ascertain(isempty(G) && isempty(W));
  ascertain(almostequal(C0,Sig));
  ascertain(almostequal(G0,Sig));
  ascertain(almostequal(W0,Sig));
  %
  disp('OK')
end  
