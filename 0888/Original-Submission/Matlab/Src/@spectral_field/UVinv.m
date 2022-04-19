
function [U,V] = UVinv(z,d) 
% UVinv - horizontal inversion  of relationship between velocity and
%         vorticity and divergence as spectral_fields.
% Usage: (U,V) = UVinv(z,d)
% Compute the velocity from the grid point values of 
% vorticity (z) and divergence (d) using the spectral
% transform and modifying the spectral coefficients.
% The grid point values are returned from a spherical 
% harmonic transform synthesis.
% But U = u *cos(lat), V = v* cos(lat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Input: 
%    z - vorticity: spectral_field object 
%      = k \cdot \nabla \times (u,v) 
%      = \frac{1}{1-\mu^2} \frac{\partial v}{\partial \lambda}
%                        - \frac{\partial u}{\partial \mu}
%    d - divergence:  spectral_field object 
%      = \nabla \cdot (u,v) 
%      = \frac{1}{1-\mu^2} \frac{\partial u}{\partial \lambda}
%                        + \frac{\partial v}{\partial \mu}
%  Output:
%    U,V - spectral_field class objects
%        In particular:
%           U.G  - field gauss_grid
%           U.gp - field grid point values
%           U.sc - field spectral coefficients
%  Local
%    xf - matrix of Fourier coefficients ordered (m,j)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by John Drake
% Method of spectral_field class: Sept 2005
%  ToDo:  add flag to see if transform coefficients are available
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ni = get(z.G,'ni');
nj = get(z.G,'nj');
mm = get(z.G,'mm');
nn = get(z.G,'nn');
kk = get(z.G,'kk');
wg = get(z.G,'wg');
yg = get(z.G,'yg');
aradius = get(z.G,'radius');
P  = get(z.G,'P');
H  = get(z.G,'H');
% Compatibility check
nnj = get(d.G,'nj');
if (nj ~= nnj) then
  error('UVinv - incompatible spectral fields on input')
else
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z =shtrana(z);   %make sure spectral coefficients are up to date
d =shtrana(d);
% Get the spectral coefficients
scz = get(z,'sc');
scd = get(d,'sc');
%
% Modify the spectral coefficients of z and d for inversion
%
%   U (\lambda_i, \mu_j )  = 
%    - \sum^{M{m=-M} \sum^{N(m)}{n=|m|,n \neq 0} \frac{a}{n(n+1)} 
%          [ i m d^m_n  P^m_n (\mu_j )} 
%           - z^m_n H^m_n (\mu_j )  ] \exp^{ i m \lambda_i}
%   V (\lambda_i, \mu_j )  = 
%    - \sum^{M{m=-M} \sum^{N(m)}{n=|m|} \frac{a}{n(n+1)} 
%          [ i m z^m_n  P^m_n (\mu_j )} 
%           + d^m_n H^m_n (\mu_j )  ] \exp^{ i m \lambda_i}
%  where H^m_n (\mu) = (1-\mu^2 ) \frac{d P^m_n (\mu)}{d \mu}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Invert U relationship
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sc = scd;  % define scratch space
scz(1,:)=0.0+0.0i;   % exclude n=0 from the sums
scd(1,:)=0.0+0.0i;   % exclude n=0 from the sums
sc(1,:) =0.0+0.0i;   % exclude n=0 from the sums
for n=1:nn
  scz(n+1,:)= scz(n+1,:) *aradius / ((n+1)*n);   
  scd(n+1,:)= scd(n+1,:) *aradius / ((n+1)*n);   
  for m=0:mm
    sc(n+1,m+1)= i * m *scd(n+1,m+1);   
  end
end
% finish transform to spectral for longitudinal derivative
xf1 =legtranOLs(sc,nj,mm,nn,kk,P);  % inverse Legendre transform 
% finish transform to spectral for latitudinal derivative
xf2 =legtranOHs(scz,nj,mm,nn,kk,H);  % inverse Legendre transform 
xf1 = -( xf1 - xf2 );  % combine Fourier coefficients for U
%
%Extend xf make xf into a hermetian array for real transform back
% Note ordering of wave numbers: 0, 1, 2, ni/2, -ni/2+1, -ni/2+2, ..,-1
% with corresponding array index:1, 2, 3, ..,ni/2+1, ni/2+2, ..........,ni
for m=1:ni/2
    xf1(ni-m+1,:) = conj(xf1(m+1,:));
end
U = spectral_field('U-velocity',z.G);
U.gp =real(ifft(xf1,ni))*ni; % inverse Fourier transform, note normalization
% But U = u *cos(lat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Invert V relationship
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scz(1,:)=complex(0.0);   % exclude n=0 from the sums
scd(1,:)=complex(0.0);   % exclude n=0 from the sums
sc(1,:) =complex(0.0);   % exclude n=0 from the sums
%  coefs of z and d already divided by n(n+1)
for m=1:mm
    n=m:nn;
    sc(n+1,m+1)= i * m *scz(n+1,m+1) ;   %(should be correct)
end
% finish transform to spectral for longitudinal derivative
xf1 =legtranOLs(sc,nj,mm,nn,kk,P);  % inverse Legendre transform 
% finish transform to spectral for latitudinal derivative
xf2 =legtranOHs(scd,nj,mm,nn,kk,H);  % inverse Legendre transform of deriv 
xf1 = -( xf1 + xf2 );  % combine Fourier coefficients for V (should be correct)
%
%Extend xf to make xf into a hermetian array for real transform back
% Note ordering of wave numbers: 0, 1, 2, ..,ni/2, -ni/2+1, -ni/2+2, ..,-1
% with corresponding array index:1, 2, 3, ..,ni/2+1, ni/2+2, ..........,ni
for m=1:ni/2
    xf1(ni-m+1,:) = conj(xf1(m+1,:));
end
V = spectral_field('V-velocity',z.G);
V.gp =real(ifft(xf1,ni))*ni;   % inverse Fourier transform, note normalization
% But V = v *cos(lat)
% all done  - U,V are new and z,d unchanged.  
% But U, V do not have active spectral coefficients
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end UVinv
