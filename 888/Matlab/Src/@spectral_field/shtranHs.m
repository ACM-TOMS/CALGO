
function f = shtranHs(f) 
% SHTRANHs Spectral H synthesis to a spectral_grid field against H.
% Compute the grid point values from the spectral coefficients  using
% a spherical harmonic H transform synthesis
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
%    H - derivatieve of associated Legendre functions ordered (nj/2,n,m)
%  Output:
%    f.sc - Output array of complex spectral coeffiecients (n,m)
%  Local
%    xf - matrix of Fourier coefficients ordered (m,j)
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
H  = get(f.G,'H');
xf =legtranOHVs(f.sc,nj,mm,nn,kk,H);      % inverse H-Legendre transform 
%make xf into a hermetian array for real transform back
for m=1:ni/2
    xf(ni-m+1,:) = conj(xf(m+1,:));
end
f.gp =real(ifft(xf,ni))*ni;           % inverse Fourier transform, note normalization

