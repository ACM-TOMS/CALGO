
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


function example_rom_nystrom_multiple_monostatic
%-----------------------------------------
% set main parameters
%-----------------------------------------

% wavenumber
kwave = 2*pi;

% setup incident plane wave
p = plane_wave(pi/2,kwave);

%-----------------------------------------
% set up scatterers
%-----------------------------------------

% Note: the code can accommodate any number of scatterers but we assume
% that they are all translations of a smaller number of possible shapes. It
% is efficient to compute the T-matrix for each shape, and reused them.

% list of possible shapes... the T-matrix will be computed for each of
% these
scatterer{1} = obstacleCircle(1);
scatterer{2} = obstacleCircle(0.5);
scatterer{3} = obstacleKite();
scatterer{4} = obstaclePinchedBall();
scatterer{5} = obstacleTrefoil();


%-----------------------------------------
% setup the solvers
%-----------------------------------------

% setup solvers with empty incident field for now
% for each type of shape
for k=1:length(scatterer)
    
    % create solver object
    solver{k} = solverNystrom(kwave,[],scatterer{k});
    
    % set Nystrom parameter
    solver{k}.setup(60)
    
end

%-----------------------------------------
% derived parameters
%-----------------------------------------

% get suggested order for the wavefunction expansion by taking the maximum
% over all possible shapes
nmax = suggestedorder(kwave,solver{1}.getRadius());
for k=2:length(scatterer)
    nmax = max(nmax,suggestedorder(kwave,solver{k}.getRadius()));
end

%-----------------------------------------
% setup the T-matrices
%-----------------------------------------

% We create the T-matrix for each possible shape

for k=1:length(scatterer)
    
    % setup T-matrix
    tmat{k} = tmatrix(nmax,kwave,solver{k},0);
    
    % print the T-matrix error estimate... based on Symmetry condition
    fprintf('T-matrix %d (%s) error estimate %0.2e\n',k,...
        class(tmat{k}.solver.scatterer),tmat{k}.error());
    
end

disp('Computed: Input/output/multiple-configuration  independent')
disp('ROM characterization of several scatterers')

% A three scatterer multiple particle configuration 
% with centers (-2, -1), (2, -1), and (0,2)
disp('Multiple config.: Choosing three scatterers (3,4,5) and their locations')


% number of scatterers
numscat = 3;

% set up arrays for scatterer types and positions... type(j) is a pointer
% into the cell array of possible shapes and sets the shape of the jth
% scatterer.
% pos(j) is the location of the jth scatterer (real part gives the x
% coordinate and imaginary part gives the y coordinate)
type = zeros(numscat,1);
pos = zeros(numscat,1);

% set scatterer types and positions
type(1) = 3;
pos(1) = -2-1i;

type(2) = 4;
pos(2) = 2-1i;

type(3) = 5;
pos(3) = 2i;

%-----------------------------------------
% solve iteratively multiple particle scattering problem
%-----------------------------------------

% create wave function expansions of plane wave at the centers of
% each of the scatterers
for j=1:numscat
    a{j} = regularwavefunctionexpansion(nmax,pos(j),p);
end

% apply the T-matrices to the incident field
for j=1:numscat
    
    % we temporarily set the origin for the T-matrix object to the position
    % of the jth scatterer... this allows the T-matrix to interact with
    % wavefunction expansions with the same origin
    tmat{type(j)}.setOrigin(pos(j));
    
    % apply the T-matrix to the incident field
    b{j} = tmat{type(j)} * a{j};
    
end

% Note: GMRES works with vectors so we need to convert our cell array of
% wavefunction expansions into a vector...
% setup an array right hand side coefficients (to pass to GMRES)
rhs = pack(b);

% set number of GMRES iterations
nitns = min(100,floor(numscat*(2*nmax+1)/2));

% Solve the linear system using GMRES
[x,flag,relres,iter,reshist] = gmres(@matrix_product,rhs,nitns,1e-8,1);

% convert coefficients into wavefunction expansions
c = unpack(x);

disp('Computed: Iteratively solved the multiple scattering model')

%-----------------------------------------
% visualize the scattered field
%-----------------------------------------

disp('Evaluating and Visualizing (Fig. 1):')
disp ('Output scattered field at 25,000 grid points (Fig. 1)')

figure(1)

% setup a grid
t=linspace(-10,10,500);
[x,y]=meshgrid(t,t);
z = x+y*1i;

% get the largest radius of the scatterers
maxrad = solver{type(1)}.getRadius();
for j=2:numscat
    maxrad = max(maxrad,solver{type(j)}.getRadius());
end

% get a mask for the scatterers
mask = abs(z-pos(1)) < 1.1 * maxrad;
for j=1:numscat
    mask = mask | abs(z-pos(j)) < 1.1 * maxrad;
end

% get the scattered field... this is just the sum of the radiating fields
% from each scatterer
scatfield = c{1}.evaluate(z,~mask);
for j=2:numscat
    scatfield = scatfield + c{j}.evaluate(z,~mask);
end

% plot the scattered field
surf(x,y,real(scatfield))
view([0 90]);
shading interp;
colorbar
title('Scattered field (in the plane [-10,10]x[-10,10]) exterior to a multiple particle configuration')

% add visualization of scatterers
hold on
pp = 2*pi*(0:99)/100;
for j=1:numscat
    [sx,sy,qx,qy]=scatterer{type(j)}.geom(pp);
    obj = plot(qx+real(pos(j)),qy+imag(pos(j)),'k-');
    set(obj,'linewidth',2);
end
hold off

%-----------------------------------------
% visualize the total field
%-----------------------------------------
disp('Evaluating and Visualizing (Fig. 2):')
disp ('Output total field at 25,000 grid points (Fig. 2)')

figure(2)

% get the total field
totalfield = scatfield + p.evaluate(z,~mask);

% plot the total field
surf(x,y,real(totalfield))
view([0 90]);
shading interp;
colorbar
title('Total field (in the plane [-10,10]x[-10,10]) exterior to a multiple particle configuration')
% add visualisation of scatterers
hold on
pp = 2*pi*(0:99)/100;
for j=1:numscat
    [sx,sy,qx,qy]=scatterer{type(j)}.geom(pp);
    obj = plot(qx+real(pos(j)),qy+imag(pos(j)),'k-');
    set(obj,'linewidth',2);
end
hold off

%-----------------------------------------
% indented function to implement the matrix
% product in GMRES
%-----------------------------------------

    function y = matrix_product(x)
        
        % convert vector of coefficients into wavefunction expansions
        c = unpack(x);
        
        % apply matrix
        for j=1:numscat
            
            % we temporarily set the origin for the T-matrix object to the position
            % of the jth scatterer... this allows the T-matrix to interact with
            % wavefunction expansions with the same origin
            tmat{type(j)}.setOrigin(pos(j));
            
            % initialize sum
            csum = regularzero(nmax,pos(j),kwave);

            % sum contributions from the other scatterers
            for i=1:numscat
                
                if i~=j
                    
                    % get the expansion of c{i} at pos{j}
                    csum = csum + regularwavefunctionexpansion(c{i},pos(j));
                    
                end
                
            end
            
            % apply the T-matrix to the sum
            d{j} = c{j} - tmat{type(j)} * csum;
            
        end
        
        % convert coefficients into a vector
        y = pack(d);
        
    end

%-----------------------------------------
% indented function to pack wavefunction
% expansion coefficients into a vector
%-----------------------------------------

    function vec = pack(a)
        
        % create an array to hold the coefficients
        vec = zeros(2*nmax+1,numscat);

        % copy the coefficients into the array
        for j=1:numscat
            vec(:,j) = a{j}.getCoefficients();
        end
        
        % reshape the array into a vector
        vec = reshape(vec,[],1);
        
    end

%-----------------------------------------
% indented function to extract wavefunction
% expansion coefficients from a vector
%-----------------------------------------

    function a = unpack(vec)
        
        % reshape the vector into an array
        vec = reshape(vec,2*nmax+1,numscat);

        % create radiating wavefunction expansions from the columns of the
        % array
        for j=1:numscat
            a{j} = radiatingwavefunctionexpansion(nmax,pos(j),kwave,vec(:,j));
        end
        
    end

end

