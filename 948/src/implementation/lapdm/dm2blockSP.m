function [cfl, dfl, sigmaPerm, p, q, rF, JNZ] = dm2blockSP(sigma, c, d, p, q)

n = size(sigma, 1);
cfl = zeros(n,1);  dfl = zeros(1,n);

sigma = sigma(p,q); c=c(p); d=d(q);

S = isfinite(sigma);
sigmatmp = sigma+1;
sigmatmp(~S) = 0;
[i,j,sij] = find(sigmatmp);
ind = sij==(d(j)+1)'-c(i);
JNZ = sparse(i(ind),j(ind),1,n,n);
[pnew,qnew,rF,~] = dmperm(JNZ);

sigmaPerm = sigma(pnew,qnew);
p = p(pnew); q = q(qnew);
JNZ = JNZ(pnew,qnew);

NumOfFineBlock = length(rF)-1;

for k=1:NumOfFineBlock
    BlockRange = rF(k):rF(k+1)-1;
    BlockSize = rF(k+1)-rF(k);
    
    sigmaBlock = sigmaPerm(BlockRange,BlockRange);
    JNZ(1:rF(k)-1,BlockRange) = -JNZ(1:rF(k)-1,BlockRange);
    
    [cfl(BlockRange),dfl(BlockRange)] = findcd(sigmaBlock,1:BlockSize);    
end

% Permute back the sparsity pattern
JNZ(p, q) = JNZ;
% Permute back fine local offsets
cfl(p) = cfl; dfl(q) = dfl;
end