
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

%-----------------------------------------
% set main parameters
%-----------------------------------------

% wavenumber
kwave = 4*pi;

% origin
center = 0;

% set scatterer shape by setting radius function and its derivative
r = @(t) 1 + 0.3 * cos(3*t);
dr = @(t) -0.9*sin(3*t);            

% set maximum for radius (used to work out T-matrix parameters)
radius = 1.3;

%-----------------------------------------
% derived parameters
%-----------------------------------------

% setup incident plane wave
p = plane_wave(0,kwave);

% compute order required
nmax = suggestedorder(kwave,radius);

%-----------------------------------------
% setup the solver
%-----------------------------------------

solver = mfsExampleSolver(kwave,[],r,dr,100,0.05);
solver.setup();

%-----------------------------------------
% setup the T-matrix
%-----------------------------------------

% setup T-matrix
tmat = tmatrix(nmax,kwave,solver,center);

% display T-matrix error check... based on Symmetry condition
fprintf('T-matrix error estimate %0.2e\n',tmat.error());
disp('Computed: Input/output independent  characterization of a sound soft trefoil')
%-----------------------------------------
% solve scattering problem
%-----------------------------------------

% create wave function expansion of plane wave
a = regularwavefunctionexpansion(nmax,center,p);

% compute wave function expansion of scattered wave using T-matrix
b = tmat * a;


disp('Computed: Output independent ROM object for the acoustic model')

disp('Evaluating and Visualizing (Fig. 1):')
disp ('Output total field at 25,000 grid points (Fig. 1)')
%-----------------------------------------
% visualize the total field
%-----------------------------------------

% setup a grid
t=linspace(-10,10,500);
[x,y]=meshgrid(t,t);
z = x+y*1i;

% get a mask for the scatterer
mask = abs(z-center) > 1.1*radius;

% plot the total field
surf(x,y,real(b.evaluate(z,mask)+p.evaluate(z,mask)));
view([0 90]);
shading interp;
colorbar
title('Total field (in the plane [-10,10]x[-10,10]) exterior to a trefoil sound soft scatterer')
