% Example for MieSolver package.
%
%   Scattering by a single scatterer with a sound soft boundary condition.
%
% Stuart C. Hawkins - 6 December 2018

% Copyright 2014-2019 Stuart C. Hawkins
% 	
% This file is part of MIESOLVER.
% 
% MIESOLVER is free software: you can redistribute it and/or modify	
% it under the terms of the GNU General Public License as published by	
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% MIESOLVER is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with MIESOLVER.  If not, see <http://www.gnu.org/licenses/>.


clear all

% setup scatterer (circle centred at 0 with unit radius and sound soft
% boundary condition).
s=scatterer(0,1,'SOFT');

% setup incident field
direction=pi/2;
kwave=5*pi;
p=plane_wave(direction,kwave);

% setup solver
pr=MieSolver(p);
pr.addScatterer(s);

% solve Mie equations
pr.solve()

% set up grid for plotting
x=linspace(-4,4,200);
[xx,yy]=meshgrid(x,x);
z=xx+1i*yy;

% plot induced field
figure(1)
val=pr.getInducedField(z);
surf(xx,yy,zeros(size(xx)),real(val))
s.show
view([0 90]);
shading interp;
axis equal
title('Induced field')
colorbar 

% plot total field
figure(2)
val=pr.getTotalField(z);
surf(xx,yy,zeros(size(xx)),real(val))
s.show
view([0 90]);
shading interp;
axis equal
title('Total field')
colorbar

% plot far field
figure(3)
tp=linspace(0,2*pi,1000);
z=exp(1i*tp);
farfield=pr.getFarfield(z);
plot(tp,real(farfield),'k-',tp,imag(farfield),'b-')
title('Far field')

% plot RCS
figure(4)
tp=linspace(0,2*pi,1000);
z=exp(1i*tp);
farfield=pr.getRcs(z);
plot(tp,farfield)
title('RCS')
