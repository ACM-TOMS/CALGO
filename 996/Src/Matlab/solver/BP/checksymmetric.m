function [symflag,Xsym,minEigVal] = checksymmetric(X,K)
    if isfield(K,'l') && K.l > 0 
        dimLP = K.l;
    else
        dimLP = 0;
    end
    minEigVal = zeros(dimLP+length(K.s),1); 
    % minEigVal = zeros(length(K.s),1); 
    [rowSize,colSize] = size(X); 
    Xsym = zeros(rowSize,colSize); 
    if dimLP > 0
        Xsym(1:dimLP) = X(1:dimLP); 
        minEigVal(1:dimLP) = X(1:dimLP); 
        count = K.l;
    else
        count = 0; 
    end        
    symflag = zeros(1,length(K.s));
    for k=1:length(K.s)
        nk = K.s(k); 
        XMat = reshape(X(count+[1:nk*nk]),nk,nk);
        symflag(k) = norm(XMat-XMat','fro') < 1e-15*(1+norm(XMat,'fro'));
        XMat = 0.5*(XMat+XMat'); 
        if nargout == 3
            minEigVal(dimLP+k) = min(eig(XMat)); 
        end
        Xsym(count+[1:nk*nk]) = reshape(XMat,1,nk*nk); 
        count = count + nk*nk; 
    end   
end
