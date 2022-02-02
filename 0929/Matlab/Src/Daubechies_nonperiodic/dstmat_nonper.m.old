function [C]=dstmat_nonper(D,N)
% [C]=dstmat_nonper(D,N)
% this function generates the quadrature matrix for the 
% Galerkin approach of the non periodic case
% Dependencies
% L_moments.m, R_moments.m, LR_partial_mom.m
%....................Grid x_1 and y_1.......................
for l=0:N-1
    x_1(l+1)=l+.5;
    y_1(l+1)=l*(N/(N-1));
end
M=(D/2);

%................PARAMETER...........
    ind=0;
    st=0;
for k=0:N-1
    
%...................................................
    if(k>=0 & k<=M-1)
        %------------------FIRST SET OF THE ROW-----------
    for m=0:M-1
    for l=0:M-1
        a(m+1,l+1)=x_1(k+l+1)^(m);
    end
    [L_mo]=L_moments(D,k);
    b(m+1)=L_mo(m+1);
end
sol=a\b' ;
        for j=0:M-1
         C(k+1,j+1)=sol(j+1);
         end
    %------------SECOND SET OF THE ROW---------------
   elseif((k>=M & k<=N-2*M+1))
        for m=0:M-1
        for l=0:M-1
        a(m+1,l+1)=l^m;%x_1(k+l+1)^(m);
        end
        [L_mo]=LR_partial_mom(D);
        b(m+1)=L_mo(m+1);
        end
        sol=a\b' ;
        for j=M+ind:M+ind+(M-1)
        C(k+1,j+1)=sol(j-(M+ind-1));
        end
         ind=ind+1;
        
        %-----------LAST SET OF THE ROW------
    else
        for m=0:M-1
        for l=0:M-1
        a(m+1,l+1)=x_1(N-M+l+1)^(m);
        end
        [R_mo]=R_moments(D,st);
        b(m+1)=R_mo(m+1);
        end
        if(st<M-1) % Due to the limitation of number of moments for boundary function
        st=st+1;
    end
    
        sol=a\b' ;
        for j=N-M:N-1     
        C(k+1,j+1)=sol(j-(N-M)+1);
        end
    end
    end
  
        
        
        

