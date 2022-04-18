function [mo]=R_moments(D,p)
% [mo]=R_moments(D,p)
% Calculate the moments of right hand boundary functions
% Dependencies
% LR_partial_mom.m, R_daubfilt.m, tmoments.m
M=D/2;  % Here M is number of vanishing moment
mom=LR_partial_mom(D);
 for k=0:M-1
     sum2=0;
     [hk_R]=R_daubfilt(k,D);
     for m=M:M+2*k
         if(p==0)
            sum2=sum2+hk_R(m+1);
        else
            tmom=tmoments(mom,-(m+1));
            sum2=sum2+hk_R(m+1)*tmom(p+1);
         end
     
     end

     for l=0:M-1
         if(k==l)
             a(k+1,l+1)=sqrt(2)-hk_R(l+1);
         else
             a(k+1,l+1)=-hk_R(l+1);
         end
     end
     b(k+1)=sum2;
 end
 sol=a\b';
 mo=sol;
