% OMEGA_LTL  LTL-factorization of Omega
%
%  [LU,LL,INFO] = OMEGA_LTL(SU,OLOW,P,Q,KO) calculates the LTL-factorization
%  Omega = L'·L of Omega which is stored in two parts, a full upper left
%  partition, SU, and a block-band lower partition, OLOW, as returned by
%  omega_build. Omega is symmetric, only the lower triangle of SU is
%  populated, and OLOW only stores diagonal and subdiagonal blocks. On exit, L =
%  [L1; L2] with L1 = [LU 0] and L2 is stored in block-band-storage in LL. INFO
%  is 0 on success. P and Q are the dimensions of the problem and KO is a vector
%  with KO(t) = number of observed values before time t.
%
%  In the complete data case KO should be 0:r:n*r. For missing values, SU and
%  OLOW are the upper left and lower partitions of Omega_o = Omega with missing
%  rows and columns removed. In this case LU and LL return L_o, the LTL-factor
%  of Omega_o.
%
%  [LU,LL,INFO,LUD,LLD] = OMEGA_LTL(SU,OLOW,P,Q,KO,SUD,OLOWD) finds in addition
%  the derivatives of LU and LL.
%
%  METHOD SKETCH: For a full block matrix A, determine k-th block row for k = n,
%  n-1,..., 1 using the partitioning of A:
%
%                |L1'  X   M' |   |L1       |   |A1  U  B'|
%                |    Lkk' Y' | · |X' Lkk   | = |U' Akk V'|
%                |         L2'|   |M   Y  L2|   |B   V  A2|
%
%  This implies that Lkk'·Lkk = Akk - Y'·Y (so Lkk is the LTL-factor of Akk -
%  Y'·Y) and also that X' = Lkk'\(U' - Y'·M). When A = Omega, LL is given by
%  this procedure, with suitable ammendments drawing on the band structure of
%  OLOW. To find LU use the partitioning of A:
%
%               |LU' M' | . |LU   | = |SU B'|
%               |    L2'|   |M  L2|   |B  A2|
%
%  which implies that LU'·LU = SU - M'·M.

function [Lu,A,info,Lud,Ad] = omega_ltl(Su,Olow,p,q,ko,Sud,Olowd)
  DIFF = nargin>5;
  LTT.LT = true; LTT.TRANSA = true;
  n = length(ko)-1;
  h = max(p,q);
  e = ko(h+1);     % order of Su
  Lu = Su;
  A = Olow;
  if DIFF
    nPar = size(Sud,3);
    Lud = Sud;
    Ad = Olowd;
  end
  for k=n:-1:h+1
    K = (ko(k)+1 : ko(k+1)) - e;
    c = ko(k-q);
    for i = min(n,k+q):-1:k
      b = ko(i-q);
      JD = K + e - c;
      JU = (b+1 : ko(k)) - c;
      J = K + e - b;
      if i>k
        I = (ko(i)+1 : ko(i+1)) - e;
        JUD = (b+1 : ko(k+1)) - c;
        JMY = 1 : ko(k+1)-b;
        A(K,JUD) = A(K,JUD) - A(I,J)'*A(I,JMY);
        if (DIFF)
          JM = 1 : ko(k)-b;
          Ad(K,JD,:) =  Ad(K,JD,:) - ata_deriv(A(I,J), Ad(I,J,:));
          Ad(K,JU,:)=Ad(K,JU,:)-atb_deriv(A(I,J),Ad(I,J,:),A(I,JM),Ad(I,JM,:));
        end
      else
        [A(K,JD), info] = ltl(A(K,JD));
        if (info~=0) return; end
        A(K,JU) = linsolve(A(K,JD), A(K,JU), LTT);
        if (DIFF)
          Ad(K,JD,:) = ltl_deriv(A(K,JD), Ad(K,JD,:));
          Ad(K,JU,:) = back_sub_deriv(A(K,JD),Ad(K,JD,:),A(K,JU),Ad(K,JU,:));
        end
      end
    end
  end
  for i = h+1:min(n,h+q)
    b = ko(i-q);
    M = b+1:e;
    I = (ko(i)+1 : ko(i+1)) - e;
    JM = 1 : ko(k)-b;
    Lu(M,M) = Lu(M,M) - A(I,JM)'*A(I,JM);
    if (DIFF) Lud(M,M,:) = Lud(M,M,:) - ata_deriv(A(I,JM), Ad(I,JM,:)); end
  end
  [Lu, info] = ltl(Lu);
  if (DIFF) Lud = ltl_deriv(Lu, Lud); end
end

function [L,info] = ltl(A)
  % Calculates LTL-factorization of A. Only uses lower triangle of A. Returns
  % info=0 on success.
  n = size(A,1);
  [L,info] = chol(A(n:-1:1,n:-1:1));
  if info~=0, L=A; return; end
  L = L(n:-1:1,n:-1:1);
end

function Ld = ltl_deriv(L, Ad)
  % Returns the derivative of the LTL-factorization of A in Ld. L should be the
  % LTL-factor of A and Ad the derivative of A.
  n = size(L,1);
  Ld = chol_deriv(L(n:-1:1, n:-1:1)', permute(Ad(n:-1:1, n:-1:1, :), [2,1,3]));
  Ld = permute(Ld(n:-1:1, n:-1:1, :), [2,1,3]);
end
