
function f = shtrans(f) 
% SHTRANA Spectral synthesis to a spectral_grid field
% Compute the grid point values from the spectral coefficients  using
% a spherical harmonic transform synthesis
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
fc =legtranOLVs(f.sc,nj,mm,nn,kk,P);      % inverse Legendre transform 
%fc =legtranMMs(f.sc,nj,mm,nn,kk,P);      % inverse Legendre transform
%make fc into a hermetian array for real transform back
for m=1:ni/2
    fc(ni-m+1,:) = conj(fc(m+1,:));
end
f.gp =real(ifft(fc,ni))*ni;           % inverse Fourier transform, note normalization

