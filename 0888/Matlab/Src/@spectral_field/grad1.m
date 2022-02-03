
function [g1,g2] = grad(f) 
% GRAD - Gradient (spherical) of a spectral_grid field.
% Usage: [g1,g2] = grad(f)
% Compute the gradient from the grid point values using the spectral
% transform and modifying the spectral coefficients.  The grid point
% values are returned from a spherical harmonic transform synthesis
% Note:  since gradient is not invariant under choice of coordinates,
%    the spherical (\theta) coordinate gradient is grad(f)/cos(\theta)
%    where \mu = sin \theta).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Input: 
%    f - spectral_field class object
%        In particular:
%           f.G  - field gauss_grid
%           f.gp - field grid point values
%           f.sc - field spectral coefficients
%  Output:
%    [g1,g2] - (vector) spectral_field object = grad(f)
%    g1 = \frac{1}{a*} \frac{\partial f}{\partial \lambda}
%    g2 = \frac{1-\mu^2}{a} \frac{\partial f}{\partial \mu}
%  Local
%    xf - matrix of Fourier coefficients ordered (m,j)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by John Drake
% Method of spectral_field class: Sept 2005
%  ToDo:  add flag to see if transform coeffecients are available
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ni = get(f.G,'ni');
nj = get(f.G,'nj');
mm = get(f.G,'mm');
nn = get(f.G,'nn');
kk = get(f.G,'kk');
wg = get(f.G,'wg');
yg = get(f.G,'yg');
P  = get(f.G,'P');
H  = get(f.G,'H');
gp = get(f,'gp');
aradius= get(f.G,'radius');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xf =fft(gp,ni)/ni;                %note normalization for Matlab FFT
% order ni= 2n real transform coming out as r0,r1,r2...i2,i1
%adjust for derivatives and multiply Fourier coefficients by the Gauss weights
% adjust for spherical coordinates
for j=1:nj
    xf1(:,j) = wg(j)*xf(:,j)/aradius;
end
xf1 = xf;   % initialization
% modify the spectral coefficients for longitudinal derivative (Eq. 5.5)
%  {\frac{1}{a*(1-mu_j^2)} \frac{\partial \xsi}{\partial \lambda} }^m_n = \sum_{j=1}^J \frac{im}{a*(1-mu_j^2)} \xsi^m (\mu_j) w_j
for m=1:mm+1
    xf1(m,:)= i * (m-1) *xf1(m,:);   % longitundal derivative
end
% finish transform to spectral for longitudinal derivative  
sc1 =legtranOLa(xf1,nj,mm,nn,kk,P);  % inverse Legendre transform 
%
%  Now work on the latitudinal derivative  (Eq. 5.8)
%  { \frac{1}{a} \frac{\partial \xsi}{ \partial \mu} }^m_n = 
%   - \sum^J_{j=1}  \xsi^m (\mu_j ) \frac{h^m_n (\mu_j )}{a*(1-\mu^2_j ) w_j
%  where H^m_n (\mu) = (1-\mu^2 ) \frac{d P^m_n (\mu)}{d \mu}
%
for j=1:nj
    xf1(:,j) = -wg(j)*xf(:,j)/aradius;
end
% finish transform to spectral for longitudinal derivative
sc2 =legtranOHa(xf1,nj,mm,nn,kk,H);  % Legendre - H analysis transform 
%
% create two new spectral_fields
%
g1 = spectral_field('grad1',f.G);
g2 = spectral_field('grad2',f.G);
g1 = set(g1,'sc',sc1);% set the new coefficients
g2 = set(g2,'sc',sc2);% set the new coefficients
g1 = shtrans(g1);    % inverse transform to get new grid point values 
g2 = shtrans(g2);    % inverse transform to get new grid point values 
% all done  - f is unchanged, but two new spectral_fields created.
% end grad
