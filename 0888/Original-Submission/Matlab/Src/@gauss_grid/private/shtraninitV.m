function [xg,yg,wg,P,H] = shtraninitV(ni,nj,mm,nn,kk) 
%Spectral object (synthesis and analysis) initialization (vectorized version)
%    ni - number of longitudes
%    nj = number of Gauss latitudes
%    mm, nn,kk are the truncation parameters
%  Output:
%    xg - Longitude points equally spaced
%    yg - Gauss points in (-1.1) indexed by j
%    wg - Gauss weights
%    P - associated Legendre functions ordered (j,n,m)
%    H - derivative of associated Legendre functions ordered (j,n,m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by John Drake
% Based on Spherical harmonic transform formulation as matrix multiply
% of Ren-Cang Lee.  
% Date: Oct. 2002
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
njo2=nj/2;
ni=2*nj;               % number of longitudes
%Calculate the Gauss weights and points
[yg,wg]= grule(nj);
%Define the uniform longitude grid, 0 <= xg < 2*pi
dx=2*pi/(ni);
xg = [0:dx:2*pi-dx];
% Fill the P array with associated Legendre functions
%USING MATLAB (UNNORMALIZED) FUNCTIONS because 'norm' is incorrect.

% Preallocate and create P

P = zeros(njo2, nn+1, nn+1);
for n=0:nn
    Pn = (legendre(n,yg))';  % Note error in Matlab normalization
    m = 0:n;
    Nmn = (-1).^(m) .* sqrt((2*n+1)/2 .* factorial(n-m)./factorial(n+m));
    for j = 1:n+1
    P(:,n+1,j) = Pn(1:njo2,j) * Nmn(j);
    end
end

%  Fill the latitudinal derivative of the associated Legendre functions:
%  H^m_n (\mu) = (1-\mu^2 ) \frac{d P^m_n (\mu)}{d \mu}
%                    = (2n+1) \epsilon^m_n P^m_{n-1} (\mu) - n \mu P^m_n (\mu).
%

% Preallocate and create H. It doesn't look like further vectorization
% will speed the program up substantially.

H = zeros(njo2, nn+1, nn+1);  % initialization
epsmn = zeros(1,1,nn+1);
for n=1:nn
  m = 0:n;
  epsmn(1,1,m+1) = sqrt(((n)^2 - (m).^2)/(4*(n)^2 -1));      
  H(:,n+1,m+1) = (2.0*n+1) * repmat(epsmn(m+1), [njo2 1 1]) .* P(:,n,m+1) - n * repmat(yg(1:njo2)', [1 1 n+1]) .* P(:,n+1,m+1);
end
%end shtraninit
