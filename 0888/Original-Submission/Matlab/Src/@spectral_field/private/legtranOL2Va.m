
function s = legtranOL2a(x,nj,mm,P) 
% Compute a Legendre transform synthesis
%  Input: 
%    x - complex Fourier coeffecients ordered (m,j) 
%    nj = number of Gauss latitudes
%    mm =nn= kk are the truncation parameters (assumed triangular)
%    P - associated Legendre functions ordered (j,n,m)
% Output:
%    s - matrix of spectral coefficients ordered (n,m)
% Local:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by John Drake
% Based on Spherical harmonic transform formulation based on open loops
% Date: Aug. 2006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
njo2=nj/2;
nn = mm;
s = zeros(nn+1,mm+1);
tic
for m=0:mm
    s(m+1:2:nn+1,m+1) = transpose((x(m+1,1:njo2) + x(m+1,nj:-1:njo2+1))*P(:,m+1:2:nn+1,m+1));
    s(m+2:2:nn+1,m+1) = transpose((x(m+1,1:njo2) - x(m+1,nj:-1:njo2+1))*P(:,m+2:2:nn+1,m+1));
end
toc       %print the time spent in analysis
