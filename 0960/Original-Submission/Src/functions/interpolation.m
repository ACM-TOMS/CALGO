      % interpolation(x,y) computes the coefficients of the Bernstein
      % polynomial p satisfying the interpolation conditions p(x(i))=y(i)
      % for i=1,...,n+1, with n the degree of the polynomial
      function coef = interpolation(x,y)
          M = bd(x); % Computes the bidiagonal decomposition of the
                     % collocation matrix of the Bernstein basis of 
                     % degree n at x(1),...,x(n+1) with high relative
                     % accuracy
                     
          [m,n] = size(M);
          if m~=n
              % TODO (change error message)
              echo('The matrix must be square'); 
              exit;
          end;
          
          coef = y;          
          for i=1:n-1
              for j=n-i+1:n
                  coef(j) = coef(j) - coef(j-1)*M(j,j-n+i);
              end
          end
          
          for i=1:n
              coef(i) = coef(i)/M(i,i);
          end
          
          for i=1:n-1
              for j=n-1:-1:i
                  coef(j) = coef(j) - coef(j+1)*M(j-i+1,j+1);
              end
          end
      end