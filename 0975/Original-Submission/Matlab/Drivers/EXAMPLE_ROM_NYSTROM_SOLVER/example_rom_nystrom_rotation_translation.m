
% Copyright 2014, 2015 Stuart C. Hawkins and M. Ganesh.
% 	
% This file is part of TMATROM.
% 
% TMATROM is free software: you can redistribute it and/or modify	
% it under the terms of the GNU General Public License as published by	
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% TMATROM is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with TMATROM.  If not, see <http://www.gnu.org/licenses/>.

clear all

% Example demonstrating the translation addition theorem for a radiating
% wavefunction expansion. We get the radiating wavefunction expansion using 
% the T-matrix for a sample scatterer.

%-----------------------------------------
% set main parameters
%-----------------------------------------

% wavenumber
kwave = 4*pi;

% origin
center = 4-2i;

% set up scatterer object
scatterer = obstacleKite();

%-----------------------------------------
% setup the solver
%-----------------------------------------

% setup solver with empty incident field for now
solver = solverNystrom(kwave,[],scatterer);
solver.setup(60)

%-----------------------------------------
% derived parameters
%-----------------------------------------

% setup incident plane wave
p = plane_wave(0,kwave);

% get radius of scatterer
radius = solver.getRadius();

% get suggested order for the wavefunction expansion
nmax = suggestedorder(kwave,radius);

%-----------------------------------------
% setup the T-matrix
%-----------------------------------------

% setup T-matrix
tmat = tmatrix(nmax,kwave,solver,0);

% display T-matrix error check... based on Symmetry condition
fprintf('T-matrix error estimate %0.2e\n',tmat.error());

disp('Computed: Characterization of a kite shaped scatterer')
%-----------------------------------------
% solve scattering problem
%-----------------------------------------

% create wave function expansion of plane wave
a = regularwavefunctionexpansion(nmax,0,p);

% compute wave function expansion of scattered wave using T-matrix
b = tmat * a;

disp('Computed: Scattering object of the kite for a fixed incident wave')

%-----------------------------------------
% verify the addition theorem by visualising....
% set up a mesh for the visualisation
%-----------------------------------------

% setup a grid
s=linspace(-3,6,300);
t=linspace(-6,3,300);
[x,y]=meshgrid(s,t);
z = x+y*1i;

%-----------------------------------------
% visualize the radiating field
%-----------------------------------------

% compute a regular wavefunction expansion of the radiating field 
% b using a new center
c = regularwavefunctionexpansion(b,center);

% get a mask for the scatterer... this is needed because the radiating
% field blows up inside the circle circumscribing the scatterer.
mask = abs(z) < 1.1*radius;

disp('Visualizing (Fig. 1): Scattering by the kite with center at zero')

% visualise the error
figure(1)
surf(x,y,zeros(size(z)),real(b.evaluate(z,~mask)))
view([0 90]);
shading interp;
colorbar
title('Scattered field (in the plane [-3,6]x[-6,3]) exterior to a sound soft kite')

% add the scatterer
hold on
solver.visualize()
hold off


% add the origins for the expansions
hold on
plot([0 real(center)],[0 imag(center)],'kx')
hold off

% store the caxis for next figure
cx = caxis;

%-----------------------------------------
% visualize the radiating field with rotated coordinates
%-----------------------------------------

% compute a regular wavefunction expansion of the radiating field 
% b using a new center
b.rotatecoordinates(pi/6);

% get a mask for the scatterer... this is needed because the radiating
% field blows up inside the circle circumscribing the scatterer.
mask = abs(z) < 1.1*radius;

disp('Visualizing (Fig. 2): Scattering by the kite in a rotated coordinate system')

% visualise the error
figure(2)
surf(x,y,zeros(size(z)),real(b.evaluate(z,~mask)))
view([0 90]);
shading interp;
colorbar
title('Rotated coordinates: scattered field (in the plane [-3,6]x[-6,3]) exterior to a sound soft kite')

% add the origins for the expansions
hold on
plot(0,0,'kx')
hold off

% set the caxis
caxis(cx)

%-----------------------------------------
% visualize the radiating field with translated origin
%-----------------------------------------

% compute a regular wavefunction expansion of the radiating field 
% b using a new center
c = regularwavefunctionexpansion(b,center);

% get a mask for the scatterer... this is needed because the radiating
% field blows up inside the circle circumscribing the scatterer.
mask = abs(z-center) < abs(center) - 1.1*radius;

disp('Visualizing (Fig. 3): Scattering by the kite at a translated center')

% visualise the error
figure
surf(x,y,zeros(size(z)),real(c.evaluate(z,mask)))
view([0 90]);
shading interp;
colorbar
title('Translated origin: scattered field (in the plane [-3,6]x[-6,3]) exterior to a sound soft kite')
% add the origins for the expansions
hold on
plot([0 real(center)],[0 imag(center)],'kx')
hold off

% set the caxis
caxis(cx)

