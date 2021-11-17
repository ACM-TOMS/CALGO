function  [alpha]=R_firstsum_alpha(m,hk_R,hk,D,M,L)
%[alpha]=R_firstsum_alpha(m,hk_R,hk,D,M,L,q);
% Function called by R_alpha.m and is calling R_alpha.m
% Dependencies
% R_alpha.m
sum1=0;
for l=0:M-1
    for k=-M+1:M
  ind1=2*m-k+1;
%ind1=2*m+k; % done by me
       if(ind1>=1 & ind1<=L-1)
        sum1=sum1+hk_R(l+1)*hk(k+(M-1)+1)*R_alpha(ind1,l,D,L);
       end

    end
end
alpha=sum1;
        
