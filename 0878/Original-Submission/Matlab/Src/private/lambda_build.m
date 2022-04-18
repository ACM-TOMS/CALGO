%LAMBDA_BUILD  Build the matrix Lambda
%
%  LAMBDA = LAMBDA_BUILD(A, r, n) builds the n·r × n·r matrix:
%
%            Lambda =   I
%                           I
%                               I
%                                   I
%                      -Ap ... -A2 -A1  I
%                          -Ap     -A2 -A1  I
%                              ...         ..  .. 
%                                  -Ap ...    -A1  I
%
%  where each Ai is r×r and A is a block row vector of the Ai, A = [A1...Ap].
%
%  [LAM, LAMD] = LAMBDA_BUILD(A, r, n, nPar) builds also the derivative of
%  Lambda.
%
%  This function is useful for testing purposes.

function [Lam, Lamd] = lambda_build(A, r, n, nPar)
  A = makecell(A);
  p = length(A);
  Lam = mat2cell(eye(r*n), repmat(r,n,1), repmat(r,n,1));
  for i=p+1:n
    for j=1:p
      Lam{i,i-j} = -A{j};
    end
  end
  Lam = cell2mat(Lam);
  if nargout > 1
    ascertain(nPar >= r*r*p); % at least the elements of A
    Lamd = zeros(r, n, r, n, r, r, p);
    for i=p+1:n
      for j=1:p
        for l=1:r
          for c=1:r
            Lamd(l,i,c,i-j,l,c,j) = -1;
          end
        end
      end
    end
    Lamd = cat(3, reshape(Lamd, r*n, r*n, r*r*p), zeros(r*n,r*n,nPar-r*r*p));
  end
end
