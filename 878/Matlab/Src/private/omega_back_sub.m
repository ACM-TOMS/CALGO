% OMEGA_BACK_SUB  Back substitution with Cholesky factor of Omega or Omega_o
%
%  X = OMEGA_BACK_SUB(Lu, Ll, Y, p, q, ko) solves L_o'·X = Y where p and q are
%  the dimensions of the ARMA model and Lo is the Cholesky factor of Omega_o,
%  stored in Lu and Ll (as returned by omega_factor). ko(t) is the number of
%  observed values before time t. 
%
%  Y may have fewer rows than L_o (which has N rows, N being the total number of
%  observations) and then L_o(1:k,1:k)·X = Y is solved, where k is the row count
%  of Y. This will save computations. There is the restriction that Y must end
%  on a block boundary (i.e. k must be equal to ko(j) for some j), and that it
%  must have at least ko(p+1) rows.
%
%  [X,Xd] = OMEGA_BACK_SUB(Lu, Ll, Y, p, q, ko, Lud, Lld, Yd) finds in addition
%  the derivative of X. Again Y (and Yd) may have less than N rows.

function [X,Xd] = omega_back_sub(Lu, Ll, Y, p, q, ko, Lud, Lld, Yd)
  DIFF = nargout>1;
  N = ko(end);
  k = size(Y,1);
  tmax = find(ko == k, 1, 'first') - 1; 
  ascertain(~isempty(tmax));
  X = Y;
  if DIFF
    nPar = size(Lld,3);
    Xd = Yd;
  end
  if isempty(X) || tmax<0, return, end
  e = size(Lu,1);
  LTT.LT = true; 
  LTT.TRANSA = true;
  h = max(p,q);
  ascertain(e == ko(h+1));
  for t = tmax:-1:h+1
    J = ko(t-q)+1:ko(t);
    K = ko(t)+1:ko(t+1);
    JL = J - ko(t-q);
    KL = K - ko(t-q);
    X(K,:) = linsolve(Ll(K-e,KL), X(K,:), LTT);
    X(J,:) = X(J,:) - Ll(K-e,JL)'*X(K,:);
    if DIFF
      Xd(K,:,:) = back_sub_deriv(Ll(K-e,KL), Lld(K-e,KL,:), X(K,:), Xd(K,:,:));
      Xd(J,:,:)=Xd(J,:,:)-atb_deriv(Ll(K-e,JL),Lld(K-e,JL,:),X(K,:),Xd(K,:,:));
    end
  end
  if e>0
    if tmax<h, e=ko(tmax+1); Lu=Lu(1:e,1:e); end
    X(1:e,:) = linsolve(Lu, X(1:e,:), LTT); 
    if DIFF
      if tmax < h, Lud=Lud(1:e,1:e,:); end
      Xd(1:e,:,:) = back_sub_deriv(Lu, Lud, X(1:e,:), Xd(1:e,:,:)); 
    end
  end
end
