function [a,b]=R_partialsum_ro(k,p,D,L,h,supp,phi)
%[a,b]=R_partialsum_ro(k,p,D,L,h,supp,phi)
% Function called by R_ro.m
% Dependencies
% R_daubfilt.m, conn.m, R_phi_origin.m, R_alpha.m
ind1=1;
M=D/2;   % Here N is number of vanishing moment
[hk_R1]=R_daubfilt(k,D);
[hk_R2]=R_daubfilt(p,D);
con_1=-conn(1,D);
    sum1=0;
    
for l=0:M-1
    for i=0:M-1
       if(l==i)
           [y,ro]=R_phi_origin(l,D,h,supp,phi);
           sum1=sum1+hk_R1(l+1)*hk_R2(i+1)*ro;
       else
           ind2=l*M+i;
      
        if(k==l & p==i)
            a(ind1,ind2)=(1-2*hk_R1(l+1)*hk_R2(i+1));
        else
            a(ind1,ind2)=-2*hk_R1(l+1)*hk_R2(i+1);
        end
    end
    end
    end
    
sum2=0;
for i=0:M-1
    for m=M:M+2*k
        sum2=sum2+hk_R2(i+1)*hk_R1(m+1)*R_alpha(m,i,D,L);
    end
end
sum3=0;
for l=0:M-1
    for q=M:M+2*p
        sum3=sum3+hk_R1(l+1)*hk_R2(q+1)*(-R_alpha(q,l,D,L));
    end
end
sum4=0;
for m=M:M+2*k
    for q=M:M+2*p
        ind=q-m+(D-1);
        if(ind>=1 & ind<=2*(D-2)+1)
        sum4=sum4+hk_R1(m+1)*hk_R2(q+1)*con_1(ind);
        end
    end
end
ro=2*(sum1+sum2+sum3+sum4);
b=ro;
