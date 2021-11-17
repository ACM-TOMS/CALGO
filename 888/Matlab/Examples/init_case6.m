
function [Ucos,Vcos,Xi,Psi,Phi,Fcor] = init_case6(nlon,nlat,X,Y);
%  Usage [Psi,Phi] = init_case6(nlon,nlat,X,Y);
%  Input:  
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
%	   Fcor   - Coriolis term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Set the Test Case 6 values:  Rossby wave
lon = X';  % note transpose so that lon corresponds to row dimension
lat = Y';
aradius = 6.37122e6; % units - m (earth radius)
Omega   = 7.292e-5;  % units 1/s (earth angular velocity)
grav    = 9.80616;   % units m/s^2 (earth gravitational acceleration)
Phi0   = grav*8000;  % normalization constant, m^2/s^2
%
%   Rossby-Haurwitz wave 4:
%   Phillips (1959), Mon. Wea. Rev., 87, 333-345
%
%     bigr = wave number
%     rnu = east-west phase velocity 
%
ir = 4;
bigr = 4.0;
bigk = 7.848d-6;
w = 7.848d-6;
rnu = (bigr*(3.0+bigr)*w - 2.0*Omega)/((1.0+bigr)*(2.0+bigr));
%
Ucos = zeros(nlon,nlat);
Vcos = zeros(nlon,nlat);
Xi = zeros(nlon,nlat);
Psi = zeros(nlon,nlat);
Phi = zeros(nlon,nlat);
Fcor = zeros(nlon,nlat);
for i=1:numel(lat)
%...Eqn. (136)
              c=cos(lat(i));
              c2=c*c;
              cr=c^ir;
              crm1=cr/c;
              crm2=crm1/c;
              c2r=cr*cr;
%
              Ucos(i) = ( aradius*w*c + aradius*bigk*crm1* (bigr*sin(lat(i))^2 -c2 )* cos(bigr*lon(i)) ) *cos(lat(i));
%...Eqn. (137)
              Vcos(i) = -aradius*bigk*bigr*crm1*sin(lat(i))*sin(bigr*lon(i))*cos(lat(i));
%...Eqn. (138)
	      Psi(i) = -aradius^2 *w*sin(lat(i))  + aradius^2 *bigk*cr*sin(lat(i))*cos(bigr*lon(i));
%...Eqn. (135)
%             Chi(i) = 0.d0
	      Xi(i) = 2.0*w*sin(lat(i)) - bigk*sin(lat(i))*cr*(bigr*bigr + 3.0*bigr+2.0)*cos(bigr*lon(i));
%...Eqn. (139)
%
	      aa = 0.5*w*(2.0*Omega + w)*c2 + 0.25*bigk^2 *c2r * ((bigr + 1)*c2 + (2.0*bigr^2 - bigr - 2.0) - 2.0* bigr^2/c2  );
%...Eqn. (141)
	      bb=2.0*(Omega+w)*bigk/((bigr+1.0)*(bigr+2.0))*cr * ((bigr^2+2.0*bigr+2.0)-(bigr+1.0)^2 *c2 );
%...Eqn. (142)
              cc=0.25*bigk^2 *c2r *((bigr+1.0)*c2 - (bigr+2.0) );
%...Eqn. (143)
              Phi(i) = Phi0 + aradius^2 * (aa + bb*cos(bigr*lon(i)) + cc*cos(2.0*bigr*lon(i)) );
%...Eqn. (140)
%              Phis(i) = 0.0;
%  Set the Coriolis term
               Fcor(i) = 2.0*Omega*sin(lat(i));
end
%  End test case 6 specification of xi, PHI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
