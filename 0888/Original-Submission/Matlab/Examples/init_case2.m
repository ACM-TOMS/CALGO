

function [Ucos,Vcos,Xi,Psi,Phi,Fcor] = init_case2(alpha,nlon,nlat,X,Y);
%  Usage [Psi,Phi] = init_case6(nlon,nlat,X,Y);
%  Input:  alpha  - rotation angle for test case 1
%          nlon   - number of longitudes
%          nlat   - number of latitudes
%          X      - meshgrid output (or slt_grid) of longitudes (lambda)
%          Y      - meshgrid output (or slt_grid) of latitudes (theta=asin(mu))
%  Output: 
%          Ucos   - u*cos velocity
%          Vcos   - v*cos velocity
%          Xi     - vorticity
%          Psi	  - stream function
%          Phi    - geopotential height (on unit sphere of radius a)  (nlon x nlat)
%	   Fcor   - Coriolis term in rotated coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Set the Test Case 2 values
lon = X';  % note transpose so that lon corresponds to row dimension
lat = Y';
aradius = 6.37122e6; % units - m (earth radius)
Omega   = 7.292e-5;  % units 1/s (earth angular velocity)
grav    = 9.80616;   % units m/s^2 (earth gravitational acceleration)
Phi0    = 2.94e4;    % normalization constant (m^2/s^2)
U0=2.0*pi*aradius/(12.0*86400.0);   %units m/s
%
Ucos = zeros(nlon,nlat); Vcos = zeros(nlon,nlat);
Xi = zeros(nlon,nlat); Psi = zeros(nlon,nlat);
Phi = zeros(nlon,nlat); Fcor = zeros(nlon,nlat);
for i=1:numel(lat)
  Ucos(i)   = U0*( cos(alpha)*cos(lat(i)) + sin(lat(i))*cos(lon(i))*sin(alpha) )*cos(lat(i));
  Vcos(i)   = -U0*sin(lon(i))*sin(alpha)*cos(lat(i));
  Psi(i)    = -aradius*U0*(sin(lat(i))*cos(alpha) - cos(lon(i))*cos(lat(i))*sin(alpha));
  Phi(i) = Phi0 - (aradius*Omega*U0 + U0*U0/2.0 )*(-cos(lon(i))*cos(lat(i))*sin(alpha) + sin(lat(i))*cos(alpha))^2;
%  Set the Coriolis term
  Fcor(i) = 2.0*Omega*( -cos(lon(i))*cos(lat(i))*sin(alpha) + sin(lat(i))*cos(alpha) );
end
%U
% pause % debug
%V
% pause % debug
%  End test case 2 specification of U,V,Xi, Psi Phi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
