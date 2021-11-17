%TEST_FIND_CGW_DERIV  Test find_CGW_deriv
%
%  TEST_FIND_CGW_DERIV prints out "OK" if the analytical derivative found with
%  find_CGW_deriv agree with numerical differentiation.

function test_find_CGW_deriv
  fprintf('TESTING FIND_CGW_DERIV... ');
  d = 0;
  for p=0:2:2
    for q=0:2:2
      for r=1:2
        d = max(d,test(p,q,r));
      end
    end
  end
  ascertain(d<1e-8);
  fprintf('max-diff=%.1e, OK\n', d)
end

function d = test(p,q,r);
  nPar = r^2*(p+q)+r*(r+1)/2;
  x = rand(nPar,1);
  d = diff_test(@fun,x,p,q,r);
end

function [f,g] = fun(x,p,q,r)
  [A,B,Sig] = theta2mat(x,p,q,r);
  if nargout==1
    [C,G,W] = find_CGW(A,B,Sig);
  else    
    [CCd,GGd,WWd] = find_CGW_deriv(A, B, Sig);
    [C,Cd] = der2array(CCd);
    [G,Gd] = der2array(GGd);
    [W,Wd] = der2array(WWd);
    g = [vec(cell2mat(Cd)); vec(cell2mat(Gd)); vec(cell2mat(Wd))];
  end
  f = [vec(cell2mat(C)); vec(cell2mat(G)); vec(cell2mat(W))];
end
