
function xf = legtranOHs(x,nj,mm,nn,kk,H) 
% Compute a Legendre H transform synthesis -- OPEN LOOP formulation
%  Input: 
%    x - complex spectral coeffecients (n,m)
%    nj = number of Gauss latitudes
%    mm, nn,kk are the truncation parameters
%    H - derviative of associated Legendre functions ordered (j,n,m)
% Output:
%    xf - matrix of Fourier coefficients ordered (m,j)
% Local:
%    S  - matrix of Fourier coeffients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by John Drake
% Based on Spherical harmonic transform formulation as an open loop 
% Date: Jan. 2003
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
njo2=nj/2;
ni= 2*nj;
%
xf=zeros(mm+1,nj);
alternatingMatrix = repmat([-1;1], ceil((nn+1)/2),1);
%tic
for m=0:mm
 xf(m+1,1:njo2) = transpose(H(:,m+1:nn+1,m+1) * x(m+1:nn+1,m+1)); %first half 
 xf(m+1,nj:-1:njo2+1) = transpose(H(:,m+1:nn+1,m+1)*(alternatingMatrix(1:nn-m+1).*x(m+1:nn+1,m+1))); %second even 
end
%toc           %print the time in matrix multiply
