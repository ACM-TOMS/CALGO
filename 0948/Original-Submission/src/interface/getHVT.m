function hvt = getHVT(sadata)
%If SWP, getHVT returns a row n-vector hvt such that (i,hvt(i)), i=1:n, 
%forms a HVT in the unpermuted signature matrix. 
%If SIP, HVT is a row n-vector of NaNs.
%hvt = getHVT(sadata)
%sadata is an object returned by daeSA.
%
%Example:
%       n = 6; G = 9.8; L = 1.0; c = 0.1;
%       sadata = daeSA(@modified2pendula,n,G,L,c);
%       hvt = getHVT(sadata)
%
%outputs
%       hvt =
%             3     2     1     6     5     4
%
%That is, the HVT of the signature matrix is in positions
%(1,3), (2,2), (3,1), (4,6), (5,5), (6,4)
%
%See also daeSA, showStruct, isSWP.
%
%Copyright 2012 Guangning Tan, Ned Nedialkov, and John Pryce
hvt  = DAESAgetHVT(sadata);
end