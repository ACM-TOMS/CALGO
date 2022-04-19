     function out = Z0n(im,in)    
           
 % This function outputs a multicomplex number with a 1 coefficient in the imaginary term imin, and a 1 in he real term. 
 % Note that for the algorithm to work m<n. Ex: Z0n(2,3) creates the multicomplex number [1,0,0,0,0,0,1,0]. 
 % It essentially outputs a zero divisor in the complex space C_n.
 
           cj=1.5;
           val=4;
           w=2;
           u=1;
           ci=0.5;

               for j=2:in
                   for i=1:im
 
                        if j-w == 1
          
                            cj=(2*cj)-1;
         
                            val=val+cj;
         
                            ci=0.5;
         
                        elseif i-u == 1
         
                            ci=2*ci;
          
                            val=val+ci;
          
                        end
      
                     w=j;
                     u=i;
        
                     if j-i == 1
                        break
                     end
         
                   end
    
               end
               
       t=zeros(1,2^in);
       t(val)=1;
       t(1)=1;
       out=multicomplex(t);   

      end
