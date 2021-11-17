      % bd(x) computes with high relative accuracy the bidiagonal 
      % decomposition of the collocation matrix of the Bernsein basis of n 
      % degree at the points x(1),...,x(n+1)     
      function M = bd(x) 
          [m,n] = size(x);
          n = n-1;
          
          M = zeros(n+1);

          
          for i=2:n+1
              auxM = (1-x(i))^n / (1-x(i-1))^(n+1);
              M(i,1) = (1-x(i-1))*auxM;
              for j=1:i-2
                  auxM = auxM*((1-x(i-1)) * (x(i)-x(i-j))) / ((1-x(i)) * (x(i-1)-x(i-j-1)));
                  M(i,j+1) = (1-x(i-j-1))*auxM;
              end
          end
          
          for j=1:n
              xj = x(j) / (1-x(j));
              for i=j+1:n+1
                  M(j,i) = xj*(n-i+2) / (i-1);
              end
          end
          
          q=1;
          M(1,1) = (1-x(1))^n;
          for i=1:n
              q = ((n-i+1) / i) * (1 / (1-x(i))) * q;
              aux = 1;
              for k=1:i
                  aux = aux*(x(i+1)-x(k));
              end
              M(i+1,i+1) = q*(1-x(i+1))^(n-i)*aux;
          end
      end       
