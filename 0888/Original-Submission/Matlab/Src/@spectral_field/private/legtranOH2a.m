
function s = legtranOH2a(x,nj,mm,H) 
% Compute a Legendre transform synthesis
%  Input: 
%    x - complex Fourier coeffecients ordered (m,j) 
%    nj = number of Gauss latitudes
%    mm = nn = kk are the truncation parameters (assumed triangular)
%    H - (real) derivative of associated Legendre functions ordered (j,n,m)
% Output:
%    s - complex matrix of spectral coefficients ordered (n,m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by John Drake
% Based on Spherical harmonic transform formulation based on open loops
% Date: Jan. 2006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
njo2=nj/2;
nn = mm;
s = zeros(nn+1,mm+1);
%tic;
for m=0:mm
%    for j=1:njo2
         j=1:njo2;
        for n=m:2:nn
          s(n+1,m+1) = s(n+1,m+1) + (x(m+1,j) - x(m+1,nj-j+1))*H(j,n+1,m+1);
        end
        for n=m+1:2:nn
          s(n+1,m+1) = s(n+1,m+1) + (x(m+1,j) + x(m+1,nj-j+1))*H(j,n+1,m+1);
        end
%    end
end
%timtota   = toc          %print the time spent in analysis
