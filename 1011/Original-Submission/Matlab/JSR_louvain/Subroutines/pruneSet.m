function [Vp, prodVp] = pruneSet(V,prodV)
% 
% [Vp, prodVp] = pruneSet(V,prodV)
%     From a set V of nv SDP matrices and a matrix prodV with nv columns,
%     generates the set Vp by deleting the majorated matrices in V. prodVp is
%     prodV from which the columns corresponding to deleted matrices have also
%     been removed.
%

nv = length(V);
nvinit = nv;
ncheck = 0;
ind1 = 1;
Vp = V;
prodVp = prodV;

while( ind1 < nv )
    ind2 = ind1+1;
    v1 = Vp{ind1};
    rem1 = false;
    
    while (ind2 <= nv)
        ind = 1:nv;
        v2 = Vp{ind2};
        
        lambda = eig(v1 - v2);
        
        if (max(lambda) <= 0) % v1 majorated by v2
            Vp = Vp(ind~=ind1);
            prodVp = prodVp(:,ind~=ind1);
            nv = nv-1;
            rem1 = true;
            break;
        elseif (min(lambda) >= 0) % v2 majorated by v1
            Vp = Vp(ind~=ind2);
            prodVp = prodVp(:,ind~=ind2);
            nv = nv-1;
        else
            ind2 = ind2+1;
        end      
    end
    
    if (~rem1)
        ind1 = ind1+1;
    end
end