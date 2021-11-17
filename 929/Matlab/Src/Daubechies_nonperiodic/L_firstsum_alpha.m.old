function  [alpha]=L_firstsum_alpha(m,hk_L,hk,D,M,L)
%[alpha]=L_firstsum_alpha(m,hk_L,hk,D,M,L)
% This function is called by L_alpha.m and is calling L_alpha.m
% Dependies
% L_alpha.m
sum1=0;
for l=0:M-1
    for k=-M+1:M
        ind1=2*m+k;
        if(ind1>=1 & ind1<=L-1)
        sum1=sum1+hk_L(l+1)*hk(k+(M-1)+1)*L_alpha(ind1,l,D,L);
    end
    end
end
alpha=sum1;
        
