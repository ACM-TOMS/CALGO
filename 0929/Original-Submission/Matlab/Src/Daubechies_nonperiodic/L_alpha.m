function [alpha]=L_alpha(m,i,D,L)
% [alpha]=L_alpha(m,i,D,L)
% Calculation of the value of left alpha functions
% Dependencies:
% L_daubfilt.m, conn.m, L_firstsum_alpha.m

M=D/2;   % Here M is number of vanishing moment
[hk_L]=L_daubfilt(i,D);
[hk]=wfilters(['db' num2str(D/2)],'r');
con_1=-conn(1,D);
  [alpha]=L_firstsum_alpha(m,hk_L,hk,D,M,L);
     sum1=alpha;
sum2=0;
for n=M:M+2*i
    for k=-M+1:M    
        ind=2*m+k-n+(D-1);
        hk(k+(M-1)+1);
        if(ind>=1 & ind<=2*(D-2)+1)          
       sum2=sum2+hk_L(n+1)*hk(k+(M-1)+1)*con_1(ind);
end
end
end
sum=2*(sum1+sum2);
alpha=sum;
