function [alpha]=R_alpha(m,i,D,L)
% [alpha]=R_alpha(m,i,D,L)
% Calculation of the value of right boundary alpha functions
% Dependencies
% R_daubfilt.m, conn.m, R_firstsum_alpha.m
M=D/2;   % Here M is number of vanishing moment
[hk_R]=R_daubfilt(i,D);
hk=wfilters(['db' num2str(D/2)],'r');
con_1=-conn(1,D);
sum1=0;
  [alpha]=R_firstsum_alpha(m,hk_R,hk,D,M,L);
     sum1=alpha;
sum2=0;
for n=M:M+2*i
    for k=-M+1:M    
      ind=n-2*m-1+k+(D-1);
   
%ind=n+2*m+1+k; %done by me
        hk(k+(M-1)+1);
        if(ind>=1 && ind<=2*(D-2)+1)          
       sum2=sum2+hk_R(n+1)*hk(k+(M-1)+1)*con_1(ind);
        end
     end
end
sum=2*(sum1+sum2);
alpha=sum;
