
function sf = legtranOLs(x,nj,mm,nn,kk,P) 
% Compute a Legendre transform synthesis -- OPEN LOOP formulation
%  Input: 
%    x - complex spectral coeffecients (n,m)
%    nj = number of Gauss latitudes
%    mm, nn,kk are the truncation parameters
%    P - associated Legendre functions ordered (j,n,m)
% Output:
%    sf - matrix of Fourier coefficients ordered (m,j)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by John Drake
% Based on Spherical harmonic transform formulation as an open loop 
% Date: Jan. 2003
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
njo2=nj/2;
ni=2*nj;
%
sf=zeros(mm+1,nj);
alternatingMatrix = repmat([1;-1], ceil((nn+1)/2),1);
j = 1:njo2;
%tic
for m=0:mm
 sf(m+1,j) = transpose(P(:,m+1:nn+1,m+1) * x(m+1:nn+1,m+1)); %first half 
 sf(m+1,nj-j+1) = transpose(P(:,m+1:nn+1,m+1)*(alternatingMatrix(1:nn-m+1).*x(m+1:nn+1,m+1))); %second even 
end
%toc       %print the time in matrix multiply
