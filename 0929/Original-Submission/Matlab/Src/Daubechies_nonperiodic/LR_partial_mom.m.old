% For calculating moment
% function called by L_moments.m and R_moments.m
%.input..............D
function [mom]=LR_partial_mom(D)
hk=wfilters(['db' num2str(D/2)],'r');
mom(1)=1;
M=(D/2);
pmax=M-1;
%........................................
for m=1:pmax
    sum3=0;
for l=0:m-1
     %..................For combination...........
        if(l~=0)
            combos=combntns(1:m,l);
            tem=size(combos,1);
        else
            tem=1;
        end
        sum1=0;
        %.....................For calculating mue...........

        for k=0:2*M-1
            sum1=sum1+hk(k+1)*(k^(m-l));
        end
        sum3=sum3+(sqrt(2)/(2*(2^(m)-1)))*tem*sum1*mom(l+1);
    end
    mom(m+1)=sum3;
end
