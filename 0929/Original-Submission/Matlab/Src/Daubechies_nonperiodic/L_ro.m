function [sol]=L_ro(D,L,h,supp,phi)
%[sol]=L_ro(D,L,h,supp,phi)
%function to calculate the vector L_ro
% Dependencies
% L_partialsum_ro.m
M=D/2;
for k=0:M-1
    for p=0:M-1
    if(k~=p)
        [a1,b1]=L_partialsum_ro(k,p,D,L,h,supp,phi);
         b(k+1)=b1;
        a(k+1,:)=a1;
    end
end
end
sol=a\b';
for k=0:M-1
    for p=0:M-1
        if(k~=p)
        ind=(k)*M+p;
            sol(ind);
        end
end
end
