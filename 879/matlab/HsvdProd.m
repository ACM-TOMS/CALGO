function B = HsvdProd(X, B, job)

%  B = HsvdProd(X, B, job)
%
%  computes the product of computes the product of an hsvdmat X
%  and a matrix B.
%
%  X       A hsvdmat.
%  B       The matrix B.
%  job     A string specifying the operation to be performed.
%
%          'ab'   B <- X*B
%          'atb'  B <- X'*B
%          'aib'  B <- X\B
%          'aitb' B <- X'\B

%  An hsvd matrix has the form
%
%      X = diag(X1, X2, ..., Xnblocks),
%
%  where
%
%      Xk = (I - uk*uk')*Sk*(I - vk*vk').
%
%  with Sk diagonal having positive diagonal entries.
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
%     sig      The diagonal elements of the Sk.
%
%  The vectors, uk, vk, and the diagonals of Sk are stacked in order in
%  the arrays X.u, X.v, and X.sig.  abs(X.bs(i)) points to the
%  beginning of the ith block and abs(X.bs(X.nblocks-1)) = n+1.  If
%  X.bs(i+1)<0, the ith block is an identity matrix.

   r1= ones(1, size(B,2));
   nblocks = X.nblocks;

   if strcmp(job,'ab')
      
      % Compute B = X*B.

      B = Hprods(X.v, B, nblocks, X.bs);
      for i=1:nblocks
         i2 = X.bs(i+1)-1;
         if i2 > 0
            i1 = abs(X.bs(i));
            B(i1:i2,:) = (X.sig(i1:i2)*r1).*B(i1:i2,:);
         end
      end
      B = Hprods(X.u, B, nblocks, X.bs);

   elseif strcmp(job,'aib')
      
      % Compute B = X\B.

      B = Hprods(X.u, B, nblocks, X.bs);
      for i=1:nblocks
         i2 = X.bs(i+1)-1;
         if i2 > 0
            i1 = abs(X.bs(i));
            B(i1:i2,:) = B(i1:i2,:)./(X.sig(i1:i2)*r1);
         end
      end
      B = Hprods(X.v, B, nblocks, X.bs);
      
   elseif strcmp(job,'atb')
      
      % Compute B = X'*B.

      B = Hprods(X.u, B, nblocks, X.bs);
      for i=1:nblocks
         i2 = X.bs(i+1)-1;
         if i2 > 0
            i1 = abs(X.bs(i));
            B(i1:i2,:) = (X.sig(i1:i2)*r1).*B(i1:i2,:);
         end
      end
      B = Hprods(X.v, B, nblocks, X.bs);
      
   elseif strcmp(job,'aitb')
      
      % Compute B = X'*\B.


      B = Hprods(X.v, B, nblocks, X.bs);
      for i=1:nblocks
         i2 = X.bs(i+1)-1;
         if i2 > 0
            i1 = abs(X.bs(i));
            B(i1:i2,:) = B(i1:i2,:)./(X.sig(i1:i2)*r1);
         end
      end
      B = Hprods(X.u, B, nblocks, X.bs);
      
   else
      error('Error in HsvdProd: Illegal Operation.');
   end
   
return

function B = Hprods(w, B, nblocks, bs)

   for i = 1:nblocks
      if bs(i+1) > 0
         i1 = abs(bs(i));
         i2 = bs(i+1)-1;
         B(i1:i2,:) = B(i1:i2,:)- w(i1:i2,1)*(w(i1:i2,1)'*B(i1:i2,:));
      end
   end
return

