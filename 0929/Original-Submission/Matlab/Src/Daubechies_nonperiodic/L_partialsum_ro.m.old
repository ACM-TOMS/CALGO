function [a,b]=L_partialsum_ro(k,p,D,L,h,supp,phi)
%[a,b]=L_partialsum_ro(k,p,D,L,h,supp,phi)
% This function is called by L_ro
% Dependies
% L_daubfilt.m, conn.m, L_phi_origin.m.
ind1=1;
M=D/2;   % Here M is number of vanishing moment
[hk_L1]=L_daubfilt(k,D);
[hk_L2]=L_daubfilt(p,D);
con_1=-conn(1,D);
sum1=0;    
for l=0:M-1
    for i=0:M-1
      
       if(l==i)
           [y,ro]=L_phi_origin(l,D,h,supp,phi);
           sum1=sum1+hk_L1(l+1)*hk_L2(i+1)*ro;
       else
             ind2=(i)*(M)+(l);
        if(k==l & p==i)
            a(ind1,ind2)=(1-2*hk_L1(l+1)*hk_L2(i+1));
        else
            a(ind1,ind2)=-2*hk_L1(l+1)*hk_L2(i+1);
        end
    end
    end
    end
    
sum2=0;
for i=0:M-1
    for m=M:M+2*k
        sum2=sum2+hk_L2(i+1)*hk_L1(m+1)*L_alpha(m,i,D,L);
    end
end
sum3=0;
for l=0:M-1
    for q=M:M+2*p
        sum3=sum3+hk_L1(l+1)*hk_L2(q+1)*(-L_alpha(q,l,D,L));
    end
end
sum4=0;
for m=M:M+2*k
    for q=M:M+2*p
        ind=m-q+(D-1);
        if(ind>=1 & ind<=2*(D-2)+1)
        sum4=sum4+hk_L1(m+1)*hk_L2(q+1)*con_1(ind);
        end
    end
end
ro=2*(sum1+sum2+sum3+sum4);
b=ro;
