%VYW_SOLVE  Solve modified vector Yule-Walker equations
%
%  S = VYW_SOLVE(A,PLU,G) returns a cell matrix S = {S0,S1,...,Sp} where Sj =
%  cov(x(t),x(t-j)). PLU is from vyw_factorize and G = {G0,G1,...,Gq} should
%  have been found by find_CGW. The Sj are obtained by solving the modified
%  vector Yule-Walker equations shown in vyw_factorize (with Gj=0 for q<j<=p in
%  case p>q).
%
%  S = VYW_SOLVE(A,PLU,Y) where Y{i} is a three dimensional array, solves, for
%  each j, the system with right hand side [Y{1}(:,:,j), ... Y{p+1}(:,:,j)] and
%  returns the result in S{1}(:,:,j),..., S{p+1}(:,:,j). This call is used when
%  calculating derivatives of S.

function S = vyw_solve(A,PLU,Y)
  A = makecell(A);
  p = length(A);
  nY = length(Y);  % nY is q+1 for original eqns and p+1 for derivative eqns
  if isstruct(Y), Y = {Y.mat}; end
  if p==0, S={Y{1}}; return, end  % NOTHING TO DO (PURE MOVING AVERAGE)
  r = size(Y{1},1);
  nrhs = size(Y{1},3);
  g = cell(p-1,nrhs);
  I = eye(r);
  % INITIALIZE g
  for i=1:p-1
    if i < nY
      g{i} = reshape(Y{i+1},r^2,nrhs);
    else
      g{i} = zeros(r^2,nrhs); 
    end
  end
  for k = 1:nrhs
    if nY > p, Y{1}(:,:,k) = Y{1}(:,:,k) + A{p}*Y{p+1}(:,:,k)'; end
    g0(:,k) = reshape(Y{1}(:,:,k)',r^2,1);
  end
  % REMOVE ROWS/COLUMNS CORRESPONDING TO UPPER-TRIANGLE FROM ROW 0/COLUMN 0
  J = [];
  for i=1:r, J = [J ((i-1)*r+i):i*r]; end
  g = [g0(J,:); cell2mat(g)];
  % SOLVE
  [perm,L,U] = deal(PLU{:});
  LT.LT = true;
  UT.UT = true;
  s = linsolve(U,linsolve(L,g(perm,:),LT),UT); % a faster s = F\g;
  % COPY s TO S0 AND FURTHER TO S{1}
  S0 = zeros(r,r,nrhs);
  j = 1;
  for i=1:r
    j1 = j+r-i;
    S0(i:r,i,:) = s(j:j1,:);
    j = j1+1;
  end
  for k=1:nrhs
    S0(:,:,k) = S0(:,:,k) + S0(:,:,k)';
  end
  S{1} = S0;
  % COPY REST OF s TO S{2},..., S{p}
  s = reshape(s(j:end,:),r,r,p-1,nrhs);
  for i=1:p-1
    S{i+1} = reshape(s(:,:,i,:),r,r,nrhs);
  end
  % FINALLY DETERMINE Sp IN S{p+1}
  if nY > p
    S{p+1} = Y{p+1};
  else
    S{p+1} = zeros(r,r,nrhs);
  end
  for k=1:nrhs
    for i=1:p
      S{p+1}(:,:,k) = S{p+1}(:,:,k) + A{i}*S{p-i+1}(:,:,k);
    end
  end
end
