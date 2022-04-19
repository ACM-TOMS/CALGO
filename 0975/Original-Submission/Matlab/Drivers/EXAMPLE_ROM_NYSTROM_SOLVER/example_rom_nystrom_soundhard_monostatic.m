
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
close all

%-----------------------------------------
% set main parameters
%-----------------------------------------

% wavenumber
kwave = 4*pi;

% origin
center = 0;

% Robin BC parameter (Sound hard case == 0)
robin_parameter = 0;

% set up scatterer object
scatterer = obstaclePinchedBall();

%-----------------------------------------
% setup the solver
%-----------------------------------------

% setup solver with empty incident field for now
solver = solverNystromRobin(kwave,[],scatterer,robin_parameter);
solver.setup(55)

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
tmat = tmatrix(nmax,kwave,solver,center);

% display T-matrix error check... based on Symmetry condition
fprintf('T-matrix error estimate %0.2e \n',tmat.error());

disp('Computed: Input/output independent  characterization of a sound-hard pinched ball')


%-----------------------------------------
% setup observation angles for monostatic
% cross section
%-----------------------------------------

disp('Solved and visualized (Fig. 3):') 
disp('1000 parameters input/output backscattering model monostatic ACS')
% set number of points
m = 1000;

% setup the grid
theta = 2*pi*(0:m-1)/m;

% get corresponding points on the unit circle
z = exp(1i*theta);

%-----------------------------------------
% solve scattering problem
%-----------------------------------------

% setup array to hold the far field
farfield = zeros(m,1);

for j = 1:m

    % setup incident plane wave with direction going 'away' from the
    % observation point
    p = plane_wave(pi+theta(j),kwave);
        
    % create wave function expansion of plane wave
    a = regularwavefunctionexpansion(nmax,center,p);
    
    % compute wave function expansion of scattered wave using T-matrix
    b = tmat * a;
    
    % compute the farfield (only in the observation direction)
    farfield(j) = b.evaluateFarField(z(j));
    
end

%-----------------------------------------
% visualize the monostatic cross section
%-----------------------------------------
figure(3)
plot(theta,10*log10(2*pi*abs(farfield).^2))
xlabel('Receiver direction angles (= - Transmitter direction angles)')
ylabel('Monostatic ACS (dB)')
