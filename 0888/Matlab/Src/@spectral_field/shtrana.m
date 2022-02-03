
function f = shtrana(f) 
% SHTRANA Spectral analysis of a spectral_grid field
% Compute the spectral coefficients from the field grid point values using
% a spherical harmonic transform analysis
%--------------------------------------------------------
%  Input: 
%    f - spectral_field class object
%        In particular:
%           f.gp - field grid point values
%           f.G - field gauss_grid
%  Imbedded data structure from the spectral_fields gauss_grid (f.G)
%    ni - number of equidistant longitudes
%    nj = number of Gauss latitudes
%    mm, nn,kk are the truncation parameters
%    yg - Gauss points in (-1.1) indexed by j
%    wg - Gauss weights
%    P - associated Legendre functions ordered (nj/2,n,m)
%  Output:
%    f.sc - Output array of complex spectral coeffiecients (n,m)
%  Local
%    fc - matrix of Fourier coefficients ordered (m,j)
%--------------------------------------------------------
% Written by John Drake
% Based on Spherical harmonic transform formulation as matrix multiply
% of Ren-Cang Lee.  
% Date: Oct. 2002
% Method of spectral_field class: Sept 2005
%--------------------------------------------------------
ni = get(f.G,'ni');
nj = get(f.G,'nj');
mm = get(f.G,'mm');
nn = get(f.G,'nn');
kk = get(f.G,'kk');
wg = get(f.G,'wg');
P  = get(f.G,'P');
gp = get(f,'gp');
fc =fft(gp,ni)/ni;                %note normalization for Matlab FFT
% order ni= 2n real transform coming out as r0,r1,r2...i2,i1
%multiply Fourier coefficients by the Gauss weights
for j=1:nj
    fc(:,j)=wg(j)*fc(:,j);
end
fc(mm+2:ni,:) = complex(0.0,0.0);
f.sc =legtranOLVa(fc,nj,mm,nn,kk,P);  % inverse (analysis) Legendre transform 
%f.sc =legtranMMa(fc,nj,mm,nn,kk,P);  % inverse (analysis) Legendre transform 


