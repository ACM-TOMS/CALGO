%%*************************************************************************
function d = blkEigBP(X,K1)
% Modified by Masakazu Kojima, 2014/01/19
    if isfield(K1,'l') && K1.l > 0
        mDim = K1.l + sum(K1.s);
        d = zeros(1,mDim); 
        d(1,1:K1.l) = X(1,1:K1.l); 
        count = K1.l; 
        eigcount = K1.l; 
    else
        mDim = sum(K1.s);
        d = zeros(1,mDim); 
        count = 0; 
        eigcount = 0; 
    end        
%   mDim = sum(K1.s); 
%   d = zeros(1,mDim); 
%   count = 0; 
%   eigcount = 0; 
    for k=1:length(K1.s)
        nk = K1.s(k);
        if (count==0) && (length(K1.s)==1) %%TKC----only 1 SDP block
           XMat = full(reshape(X,nk,nk));
        else
           XMat = full(reshape(X(count+[1:nk*nk]),nk,nk));
        end
        XMat = (XMat + XMat')/2;
        eigval = eig(XMat);
        d(eigcount+[1:nk])  = eigval;
        count = count + nk*nk;
        eigcount = eigcount + nk;
    end
end
%%*************************************************************************