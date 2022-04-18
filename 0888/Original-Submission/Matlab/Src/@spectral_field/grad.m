
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
%    g1 = \frac{1}{a} \frac{\partial f}{\partial \lambda}
%    g2 = \frac{1-\mu^2}{a} \frac{\partial f}{\partial \mu}
%  Local
%    x - matrix of spectral coefficients ordered (m,n)
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
%
% create two new spectral_fields
%
g1 = f;
g2 = f;
% Get the spectral coffiecients of the input field (analysis)
g1 = shtrana(g1);
x = get(g1,'sc');
% clean up the documentation
% modify the spectral coefficients for longitudinal derivative (Eq. 5.5)
%  {\frac{1}{a} \frac{\partial \xsi}{\partial \lambda} }^m_n 
%      = \sum_{j=1}^J \frac{im}{a} \xsi^m (\mu_j) w_j
%      = \frac{1}{a} \sum_{m,n} i m \xsi^m_n P^m_n(\mu) e^{i m \lambda}
for m=0:mm
    sc1(m+1,:)=  i * m *x(m+1,:)/aradius;   % longitundal derivative (i=sqrt(-1))
end
g1 = set(g1,'sc',sc1);
%
%  Now work on the latitudinal derivative  (Eq. 5.8)
%  { \frac{1-\mu^2}{a} \frac{\partial \xsi}{ \partial \mu} }^m_n 
%   = - \sum^J_{j=1}  \xsi^m (\mu_j ) \frac{H^m_n (\mu_j )}{(1-\mu^2_j ) w_j
%   = \frac{1}{a} \sum_{m,n} \xsi^m_n H^m_n(\mu) e^{i m \lambda}
%  where H^m_n (\mu) = (1-\mu^2 ) \frac{d P^m_n (\mu)}{d \mu}
%
    sc2 = x/aradius;
% finish transform to spectral for longitudinal derivative
g2 = set(g2,'sc',sc2);% set the new coefficients
g1 = shtrans(g1);    % inverse transform to get new grid point values 
g2 = shtranHs(g2);   % inverse H- transform to get new grid point values 
%g2 = shtrana(g2);    % fill in the spectral coefficients of the gradient
% all done  - f is unchanged, but two new spectral_fields created.
% end grad
