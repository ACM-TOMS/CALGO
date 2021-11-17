
function [U,V,F] = init_case1(alpha,nlon,nlat,X,Y);
%  Usage [U,V,F] = init_case1(alpha,nlon,nlat,X,Y);
%  Input:  alpha  - rotation angle for test case 1
%          nlon   - number of longitudes
%          nlat   - number of latitudes
%          X      - meshgrid output (or slt_grid) of longitudes (lambda)
%          Y      - meshgrid output (or slt_grid) of latitudes (theta=asin(mu))
%  Output: U	  - velocity in in lambda direction  (nlon x nlat) = u*cos(lat)
%          V      - velocity in the theta direction  (nlon x nlat) = v*cos(lat)
%          F      - geopotential height (on unit sphere?)  (nlon x nlat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Set the Test Case 1 values
lon = X';  % note transpose so that lon corresponds to row dimension
lat = Y';
% alpha
thetac = 0.0; lamc = 3.0*pi/2.0; bigR = 1/3;
U0=2.0*pi/12.0; H0=0.1; 
%
U = zeros(nlon,nlat); V = zeros(nlon,nlat);
for i=1:numel(lat)
  U(i) = cos(lat(i))*U0*( cos(alpha)*cos(lat(i)) + sin(lat(i))*cos(lon(i))*sin(alpha) );
  V(i) = -cos(lat(i))*U0*sin(lon(i))*sin(alpha);
end
%U
% pause % debug
%V
% pause % debug
%  Define an initial test field to be advected on grid interior
% Careful about the branch of the arcsine and alpha = pi/2
if (alpha == 0.0)
  rlata = thetac;
  rlona = lamc;
else
  tst = sin(thetac)*cos(-alpha) - cos(thetac)*cos(lamc)*sin(-alpha);
  if (tst > 1.0) 
    rlata = pi/2.0;
  elseif (tst <-1.0)
    rlata = -pi/2.0;
  else
    rlata = asin(tst);
  end
  tst = cos(rlata);
  if (tst == 0.0) 
    rlona = 0.0;
  else
    tst = sin(lamc)*cos(thetac)/tst;
    if (tst > 1.0)
      rlona =pi/2.0;
    elseif (tst < -1.0) 
      rlona = -pi/2.0;
    else
      rlona = asin(tst);
    end
end
tst = cos(-alpha)*cos(lamc)*cos(thetac) + sin(-alpha)*sin(lamc);
if (tst < 0.0)
  rlona = pi - rlona;
end
end
%
F = zeros(nlon,nlat);
for i=1:numel(lat)
  r = acos( sin(rlata)*sin(lat(i)) + cos(rlata)*cos(lat(i))*cos(lon(i)-rlona) );
  if ( r<bigR )
    F(i) = (H0/2.0 )*(1.0 + cos(pi*r/bigR));
  end
end
%  End test case 1 specification of U,V, PHI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
