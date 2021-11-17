
function g = curl(U,V) 
% CURL - horizontal curl  (spherical) of a vector spectral_field.
% Usage: g = curl(U,V)
% Compute the curl from the grid point values using the spectral
% transform and modifying the spectral coefficients.  The grid point
% values are returned from a spherical harmonic transform synthesis
%--------------------------------------------------------
%  Input: 
%    U,V - spectral_field class objects
%        In particular:
%           U.G  - field gauss_grid
%           U.gp - field grid point values
%           U.sc - field spectral coefficients
%  Output:
%    g - spectral_field object 
%      = k \cdot \nabla \times (U,V) 
%      = \frac{1}{a(1-\mu^2)} \frac{\partial V}{\partial \lambda}
%        - \frac{1}{a} \frac{\partial U}{\partial \mu}
%  Local
%    xf - matrix of Fourier coefficients ordered (m,j)
%--------------------------------------------------------
% Written by John Drake
% Method of spectral_field class: Sept 2005
%  ToDo:  add flag to see if transform coeffecients are available
%--------------------------------------------------------
ni = get(U.G,'ni');
nj = get(U.G,'nj');
mm = get(U.G,'mm');
nn = get(U.G,'nn');
kk = get(U.G,'kk');
wg = get(U.G,'wg');
yg = get(U.G,'yg');
P  = get(U.G,'P');
H  = get(U.G,'H');
aradius = get(U.G,'radius');
gpu = get(U,'gp');
gpv = get(V,'gp');
%--------------------------------------------------------
xfu =fft(gpu,ni)/ni;                %note normalization for Matlab FFT
xfv =fft(gpv,ni)/ni;                %note normalization for Matlab FFT
% order ni= 2n real transform coming out as r0,r1,r2...i2,i1
%adjust for derivatives and multiply Fourier coefficients by the Gauss weights
for j=1:nj
    xfu(:,j)=wg(j)*xfu(:,j)/(aradius*(1.0-yg(j)^2));
    xfv(:,j)=wg(j)*xfv(:,j)/(aradius*(1.0-yg(j)^2));
end
%
% Modify the Fourier coefficients of U and V for curl
%
%  { curl (U,V) }^m_n = 
%    - \sum^J_{j=1} [ i m V^m ( \mu_j ) P^m_n (\mu_j )} 
%                 + U^m ( \mu_j ) H^m_n (\mu_j )  ] \frac{w_j}{a (1-\mu^2_j )}
%  where H^m_n (\mu) = (1-\mu^2 ) \frac{d P^m_n (\mu)}{d \mu}
%
for m=0:mm
    xfv(m+1,:)= i * m *xfv(m+1,:);   % longitundal derivative (i=sqrt(-1))
end
% finish transform to spectral for longitudinal derivative
sc1 =legtranOLa(xfv,nj,mm,nn,kk,P);  % inverse Legendre transform 
% finish transform to spectral for latitudinal derivative
sc2 =legtranOHa(xfu,nj,mm,nn,kk,H);  % inverse Legendre transform 
%
% create new spectral_field with coeffients as the sum of the parts
%
%sc1 = -(sc1 + sc2);  % is the sign screwed up?
sc1 = (sc1 + sc2);  % works with SWE
g = spectral_field('curl',U.G);
g = set(g,'sc',sc1);% set the new coefficients
g = shtrans(g);    % inverse transform to get new grid point values 
% all done  - U,V are unchanged, but a new spectral_field created.
% end curl
