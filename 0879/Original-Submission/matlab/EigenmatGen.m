function A = EigenmatGen(n, nblocks, yident, zident)
   
%  A = EigenmatGen(n, nblocks, yident, zident)
%
%  gaenerates an eigenmat.
%
%  n       The order of the eigenmat.
%  nblocks The number of blocks in Z
%  yident  If nonzero, the hsvdmat Y is an identity.
%  zident  If nonzero, the hsvdmat Z is an identity.
%
%  EigenmatAlloc allocates zeroed memory for an eigenmat A.
%  In addition it initializes A.n, A.Z.n, A.Z.nblocks, A.Z.bs(0)
%  A.Z.bs(n), A.y.n, A.Y.nblocks, and A.Y.bs.

%  An eigenmat has the form
%
%      A = Y*Z*(L - shift*I)*Z^-1*Y^-1.
%
%  Here L is a block triangular with 1x1 blocks containing the
%  real eigenvalues of A, 2x2 blocks conaining the complex
%  eignevalues of A in the form
%
%               |  mu  nu |
%               |         |,
%               | -nu  mu |
%
% and Jordon blocks of order k having the form (for k=3)
%
%               | lambda     eta_1      0      |
%               |                              |
%               |   0       lambda     eta_2   |.
%               |                              |
%               |   0         0       lambda   |
%
%  Y and Z are hsvd transformations, with Y having only one block
%  (see below).
%
%  An eigenmat has the following fields.
%
%     n        The order of the matrix.
%     eig(:)   Contain the  eigenvalues of A and their
%     type(:)  types according to the following scheme.
%              Real eigenvalue
%                 type(i)=1      eig(i)=the eigenvalue.
%              Complex conjugate eigenvalues
%                 type(i)=2      eig(i)=mu
%                 type(i)=3      eig(i)=nu
%              Jordan block of order k
%                 type(i)=-k     eig(i)=lambda
%                 type(i+j)=-1   eig(i+j)=eta_j
%     Y, Z     The outer and inner hsvd  transformations.
%
%  An hsvd matrix has the form
%
%      X = diag(X1, X2, ..., Xnblocks),
%
%  where
%
%      Xk = (I - uk*uk^T)*Sk*(I - vk*vk^t).
%
%  with Sk diagonal having positive diagonal entries.  The vectors,
%  uk, vk, and the diagonals of Sk are stacked in order in the arrays
%  X%u, X%v, and X%sig.  abs(X%bs(i)) points to the beginning of the ith
%  block and abx(X.bs(X.nblocks-1)) = n.  If X.bs(i+1)<0, the ith block
%  is an identity matrix.
%  
%  An hsvdmat has the following fields.
%
%     n        The order of the matrix.
%     nblocks  The number of blocks in the hsvdmat
%     bs       abs(bs(i)) is the index of the start of the i-th block.
%              abs(bs(nblocks+1)) = n+1.  If bs(i+1)<0, the i-th block
%              is an identity.
%     u        The vector generating the left Householder transformation.
%     v        The vector generating the right Householder transformation.
%     sig      The singular values.

   if yident
      Y = struct('n', n, 'nblocks', 1, 'bs', [1;-(n+1)], ...
                 'sig', [], 'u', [], 'v', []);
   else
      Y = struct('n', n, 'nblocks', 1, 'bs', [1; n+1], ...
                 'sig', zeros(n,1), 'u', zeros(n,1), 'v', zeros(n,1));
   end
   
   if zident
      Z = struct('n', n, 'nblocks', 1, 'bs', [1;-(n+1)], ...
                 'sig', [], 'u', [], 'v', []);
   else
      Z = struct('n', n, 'nblocks', nblocks, ...
                 'bs', [1;zeros(nblocks-1,1);(n+1)], ...
                 'sig', zeros(n,1), 'u', zeros(n,1), 'v', zeros(n,1));
   end
   
   A = struct('n', n, 'eig', zeros(n,1), 'type', zeros(n,1), 'Z', Z, 'Y', Y);
   
return
