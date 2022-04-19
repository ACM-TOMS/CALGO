function [eig, x, y, cond] = EigenmatVecs(A, eignum, job)
%
%  [eig, x, y, cond] = EigenmatVecs(A, eignum, job)
%
%  computes an eigenvalue the corresponding left or 
%  right eigenvectors as specified job.  If both
%  are computed,  EigenmatVecs also returns the condition
%  number of the eigenvalue.  The eigenvectors are scaled
%  to have Euclidean norm 1.
%
%  A        The eigenmat whose vectors are to be computed.
%  eignum   The position in A.eig of the eigenvalue.
%  job      A string specifying what to compute.
%
%           "r"  Compute the right eigenvector.
%           "l"  Compute the left eigenvector.
%           "b"  Compute both and the condition number.
%
%                (Note: For Jordan blocks, principal vectors
%                are computed and -1 is returned for the
%                condition number.)
%
%  eig      The eigenvalue.
%  x(:)     The right eigenvector.
%  y(:)     The left eigenvector.
%  cond     The condition number of the eigenvalue.
%           (or -1, if the eigenvalue belongs to a
%            Jordan block)

   imag1 = sqrt(-1);

   x = zeros(A.n,1);
   y = zeros(A.n,1);
   cond = 0;
   
   if (eignum > A.n)
      error('Error in EigenmatVecs: eignum > n.')
   end
   
   if (~strcmp(job, 'b') && ~strcmp(job, 'r') && ~strcmp(job, 'l'))
      error('Error in EigenmatVecs: Illegal operator.')
   end
   
   if (A.type(eignum) <=1) 

      % Real egenvalue.

      if (A.type(eignum) == 1 || A.type(eignum) < -1) 
          eig = A.eig(eignum);
      elseif (A.type(eignum) == -1)
          found = 0;
          for i=eignum-1:-1:1
              if (A.type(i) < -1)
                  eig = A.eig(i);
                  found = 1;
                  break;
              elseif (A.type(i) ~= -1) 
                 error('Error in EigenmatVecs: Illegal type in Jordan block.')
              end
          end
          if (found == 0)
              error('Error in EigenmatVecs: Illegal type in Jordan block.')
          end
      else
          error('Error in EigenmatVecs: Illegal type.')
      end
          
      if (strcmp(job, 'r') || strcmp(job, 'b'))

          % Compute the right eigenvector.

          x(eignum) = 1;
          x = HsvdProd(A.Z, x, 'ab');
          x = HsvdProd(A.Y, x, 'ab');
          if (A.type(eignum) == 1)
              x = x/norm(x);
          end
      end
      if (strcmp(job, 'l') || strcmp(job, 'b'))

          % Compute the left eigenvector.

          y(eignum) = 1;
          y = HsvdProd(A.Z, y, 'aitb');
          y = HsvdProd(A.Y, y, 'aitb');
          if (A.type(eignum) == 1)
              y = y/norm(y);
          end
      end

      if (strcmp(job, 'b'))

         % Compute the condition number.

         if (A.type(eignum) == 1)
            cond = 1/abs(y'*x);
         else
            % a Jordan block    
            cond = -1; 
         end
      end

   elseif (eignum == A.n)
      
      error('Error in EigenmatVecs: eignum is too large.')
      
   elseif (A.type(eignum) == 2)
      
      if (A.type(eignum+1) ~= 3) 
         error('Error in EigenmatVecs: type 2 should be followed by 3.')
      end
      
      % Compute the eigenvalue.

      eig = complex(A.eig(eignum), A.eig(eignum+1));
      
      if (strcmp(job, 'r') || strcmp(job, 'b'))
  
         % Compute the right eigenvector.

         x(eignum) = 1;
         x(eignum+1) = imag1;
         x = HsvdProd(A.Z, x, 'ab');
         x = HsvdProd(A.Y, x, 'ab');
         x = x/norm(x);
      end
         
      if (strcmp(job, 'l') || strcmp(job, 'b'))

         % Compute the left eigenvector.

         y(eignum) = 1;
         y(eignum+1) = imag1;
         y = HsvdProd(A.Z, y, 'aitb');
         y = HsvdProd(A.Y, y, 'aitb');
         y = y/norm(y);
      end
         
      if (strcmp(job, 'b'))

         % Compute the condition number.

         cond = 1/abs(y'*x);
      end
         
   else
      A.type(eignum:eignum+1)

      error('Error in EigenmatVecs: Illegal type.')
   end
return
