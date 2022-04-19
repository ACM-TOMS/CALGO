      function [r,k] = powerDC(x,n)
      
        if n==0
            r = 1;
        else
            i = n-1;
            r = x;
            y=x;
            k = 0;

            while i>0
                if rem(i,2)==1 
                    r = r.*y;
                    k = k+1;                 
                end
                i = floor(i/2);
                if i~=0
                    y=y.*y;
                    k = k+1;
                end
            end % while
      end % powerDC