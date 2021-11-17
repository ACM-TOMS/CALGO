%VYW_FACTORIZE  Set up and LU-factorize modified vector Yule-Walker equations
%
%  PLU = VYW_FACTORIZE(A) returns the LU-factorization of a matrix F of order
%  r^2·p-r(r-1)/2 which may be used to solve the modified Yule-Walker equations:
%
%     S0 - A1*S1' + A2*S2' + A3*S3' + ... + Ap*Sp'     = G0   (1.0)
%     S1 - A1*S0  + A2*S1' + A3*S2' + ... + Ap*S(p-1)' = G1   (1.1)
%     S2 - A1*S1  + A2*S0  + A3*S1' + ... + Ap*S(p-2)' = G2   (1.2)          (1)
%     S3 - A1*S2  + A2*S1  + A3*S0  + ... + Ap*S(p-3)' = G3   (1.3)
%     ...
%     Sp - A1*S(p-1) + A2*S(p-2) + ........ + Ap*S0    = Gp   (1.p)
%
%  where Gj = cov(y(t),x(t-j)) and Sj = cov(x(t),x(t-j)). The system is obtained
%  by substituting Sp given by (1.p) into (1.0) and then making modifications to
%  take into account that S0 is symmetric. The result is a linear system
%
%                          F·s = g
%
%  where s and g are column vectors, g obtained from the Gj and s containing the
%  elements of S0, S1,..., S(p-1) in column order with upper triangle elements
%  of S0 removed (i.e. s' = [S0(1,1),..., S0(1,r), S0(2,2),..., S0(r,r), 
%  S1(1,1),..., S1(1,r), S1(2,1),..., S(p-1)(r,r)])
%
%  PLU returns a cell array {perm, L, U} where perm specifies a permutation
%  matrix P such that P·L·U = F. On entry A should be [A1...Ap].

function PLU = vyw_factorize(A)
  A = makecell(A);
  p = length(A);
  if p==0, PLU={}; return, end  % NOTHING TO DO (PURE MOVING AVERAGE)
  r = size(A{1},1);
  F0r = cell(1,p-1);
  F0c = cell(p-1,1);
  F = cell(p-1,p-1);
  I = eye(r);
  % SET F TO IDENTITY MATRICES, F0c TO ZERO AND INITIALIZE F0r AND F00
  for i=1:p-1
    for j=1:p-1, F{i,j} = zeros(r^2,r^2); end
    F{i,i} = eye(r^2);
    F0c{i} = zeros(r^2,r^2);
    F0r{i} = -kron(A{i},I) - kron(A{p},A{p-i});
  end
  F00 = eye(r^2) - kron(A{p},A{p});
  % MAIN LOOP. LOOP OVER SUBMATRICES (K,J), BLOCK-ROWS (i) AND BLOCK-COLUMNS (j)
  K = 1:r;
  J = 1:r:r^2;
  for k = 1:r
    for i=1:p-1
      for j=1:i-1   % SUBTRACT AHAT MATRICES FROM F
        F{i,i-j}(K,K) = F{i,i-j}(K,K) - A{j};
      end
      for j=i+1:p   % SUBTRACT ATILDA MATRICES FROM F
        F{i,j-i}(K,J) = F{i,j-i}(K,J) - A{j};
      end
      F0c{i}(K,K) = F0c{i}(K,K) - A{i}; % SUBTRACT AHATi FROM COLUMN 0
      F0c{i}(K,J) = F0c{i}(K,J) - A{i}; % SUBTRACT ATILDAi FROM COLUMN 0
    end
    % SUBTRACT ATILDAp·AHATp FROM F00 & ADD TO DIAGONAL OF F00
    K1 = 1:r;
    for k1 = 1:r
      F00(K,K1) = F00(K,K1) - A{p}(:,k1)*A{p}(k,:);
      K1 = K1+r;
    end
    F00(K(k),K(k)) = F00(K(k),K(k)) + 1;
    K = K+r;
    J = J+1;
  end
  % CHANGE F0r, F0c, AND F FROM CELL-ARRAYS TO ORDINARY MATRICES
  F0r = cell2mat(F0r);
  F0c = cell2mat(F0c);
  if p <= 1, F0r = zeros(r^2,0); F0c = zeros(0,r^2); end
  F = cell2mat(F);
  % REMOVE ROWS/COLUMNS CORRESPONDING TO UPPER-TRIANGLE FROM ROW 0/COLUMN 0
  J = [];
  for k=1:r, J = [J ((k-1)*r+k):k*r]; end
  F = [
    F00(J,J) F0r(J,:)
    F0c(:,J)    F    ];
  [L,U,P] = lu(F); [perm,dum] = find(P');
  PLU = {perm,L,U};
  if ~all(diag(U)), perm(1)=0; end
end
