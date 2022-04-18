
function x = legtranOL2s(s,nj,mm,P) 
% Compute a Legendre transform synthesis -- OPEN LOOP formulation
%  Input: 
%    s - complex spectral coeffecients (n,m)
%    nj = number of Gauss latitudes
%    mm = nn = kk are the truncation parameters (assumed triangular)
%    P - associated Legendre functions ordered (j,n,m)
% Output:
%    x - complex matrix of Fourier coefficients ordered (m,j)
% Local:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by John Drake
% Based on Spherical harmonic transform formulation as an open loop 
% Date: Aug. 2006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
njo2=nj/2;
nn = mm;
ni=2*nj;
%
x=zeros(mm+1,nj);
%tic;
for m=0:mm
    for j=1:njo2
        for n = m:nn
           x(m+1,j) = x(m+1,j) +  s(n+1,m+1)*P(j,n+1,m+1); %first half
        end
        for n=m:2:nn
           x(m+1,nj-j+1) = x(m+1,nj-j+1) +  s(n+1,m+1)*P(j,n+1,m+1); %second even
        end
        for n=m+1:2:nn
           x(m+1,nj-j+1) = x(m+1,nj-j+1) -  s(n+1,m+1)*P(j,n+1,m+1); %second odd
        end
    end
end
%timtots = toc            %print the time in matrix multiply
