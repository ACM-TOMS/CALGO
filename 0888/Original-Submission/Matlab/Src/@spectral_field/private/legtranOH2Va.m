
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
%tic
for m=0:mm
    s(m+1:2:nn+1,m+1) = transpose((x(m+1,1:njo2) - x(m+1, nj:-1:njo2+1))*H(:,m+1:2:nn+1,m+1));
    s(m+2:2:nn+1,m+1) = transpose((x(m+1,1:njo2) + x(m+1, nj:-1:njo2+1))*H(:,m+2:2:nn+1,m+1));
end
%toc         %print the time spent in analysis
