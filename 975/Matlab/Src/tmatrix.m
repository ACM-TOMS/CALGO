% T-matrix of a scatterer
%
%   T = tmatrix(n,k,s,x) returns the T-matrix of order n, wavenumber k, and
%   origin x. The scatterer details are contained in the object S which
%   must be an instance of the 'solver' class.
%
% Also:
%
%   val = T.error() gives a measure of the error in the T-matrix based on a
%   symmetry relation. See Ganesh and Hawkins ANZIAM J. Vol 51 
%   Pages C215--C230 (2010) for details.
% 
%   T.setOrigin(x) sets the origin of the T-matrix to x. This is a virtual
%   origin that is only used in interactions with wave functions. This does
%   not change the T-matrix.
%
% Example:
%
%   p = plane_wave(0,k);
%   u = regularwavefunctionexpansion(n,0,p);
%   T = tmatrix(n,k,s,0);
%   v = T * u;
%
%   Now v is a radiating wavefunction expansion for the scattered field
%   induced by the plane wave p.
%
% See also: regularwavefunctionexpansion, radiatingwavefunctionexpansion,
% plane_wave, point_source.
%
% Stuart C. Hawkins - 12 November 2014

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


classdef tmatrix < handle
    
    properties
        
        % given properties
        order
        kwave
        solver
        origin
        
        % derived properties
        matrix
        
    end
    
    methods

        %-----------------------------------------
        % constructor
        %-----------------------------------------

        function self = tmatrix(order,kwave,solver,origin)
            
            % set default for origin
            if nargin < 4
                origin = 0;
            end
            
            % check that solver is of the correct type
            if ~isa(solver,'solver')
                error('solver must be an instance of the solver class')
            end
            
            % - - - - - - - - - - - - - - - - - 
            % set properties
            % - - - - - - - - - - - - - - - - - 
            
            % set given properties
            self.order = order;
            self.kwave = kwave;
            self.solver = solver;
            self.origin = origin;
            
            % - - - - - - - - - - - - - - - - - 
            % setup the quadrature points
            % - - - - - - - - - - - - - - - - - 

            % get quadrature points and weights
            points = pi*(0:2*self.order+1)/(self.order+1);
            weight = pi/(self.order+1);
            
            % ensure these are column vector
            points = points(:);
            
            % - - - - - - - - - - - - - - - - - 
            % solve the scattering problems
            % - - - - - - - - - - - - - - - - - 

            % setup a cell array of incident fields
            for n = -self.order:self.order
                inc{n + self.order + 1} = regularwavefunction2d(n, ...
                    self.kwave,self.origin);
            end
            
            % set the incident field in the solver
            self.solver.setIncidentField(inc);
            
            % solve the scattering problems
            self.solver.solve();
            
            % get the farfield
            farfield = self.solver.getFarField(points,1:2*self.order+1);
            
            % assume the farfield was computed with origin 0... we need
            % to adjust for the modified origin
            if self.origin ~= 0
               
                % Use (21) in Dufva et al, Progress in Electromagnetics
                % Research B, Vol 4, 79-99, 2008
                
                % compute scaling factor
                sigma = exp(1i*self.kwave*real(conj(self.origin)*exp(1i*points)));
                
                % transform the far field
                farfield = spdiags(sigma(:),0,length(sigma),length(sigma)) * farfield;
                
            end
            
            % - - - - - - - - - - - - - - - - - 
            % compute the T-matrix
            % - - - - - - - - - - - - - - - - - 

            % get the order of the wavefunctions in a row vector
            n = (-self.order:self.order);

            % compute the T-matrix using (16) in Ganesh and Hawkins
            % ANZIAM J. 51 C215--C230 (2010)
            self.matrix = weight * self.scaling(n) * exp(1i*points*n)' * farfield;
            
        end
            
        %-----------------------------------------
        % scaling function
        %-----------------------------------------

        % gives the coefficient in (16) Ganesh and Hawkins
        % ANZIAM J. 51 C215--C230 (2010)
        
        function val = scaling(self,n)
           
            vec = 0.25*(1+1i)*sqrt(self.kwave/pi)*1i.^abs(n);
            
            val = spdiags(vec(:),0,length(vec),length(vec));
            
        end            
        
        %-----------------------------------------
        % error check
        %-----------------------------------------

        % Check the symmetry relation (17) in Ganesh and Hawkins
        % ANZIAM J. 51 C215--C230 (2010)

        function val = error(self)
           
            val = max(max(abs(self.matrix + self.matrix' ...
                + 2 * self.matrix' * self.matrix)));
            
        end
                    
        %-----------------------------------------
        % multiply tmatrix x regular wave expansion
        %-----------------------------------------

        function val = mtimes(self,expansion)
            
            % - - - - - - - - - - - - - - - - - 
            % check the T-matrix and the wave 
            % expansion are compatible
            % - - - - - - - - - - - - - - - - - 

            if ~isa(expansion,'regularwavefunctionexpansion')
                
                error('expansion must be a regularwavefunctionexpansion')
                
            end
            
            if self.kwave ~= expansion.kwave
                
                error('T-matrix and expansion wavenumbers do not match.')
                
            end
            
            if self.origin ~= expansion.origin
                
                error('T-matrix and expansion centers do not match.')
                
            end
           
            if self.order ~= expansion.order
                
                error('T-matrix and expansion orders do not match.')
                
            end
            
            % - - - - - - - - - - - - - - - - - 
            % do product 
            % - - - - - - - - - - - - - - - - - 

            % create a radiating wave function expansion with coefficients
            % obtained by matrix multiplication with the T-matrix
            val = radiatingwavefunctionexpansion(self.order,self.origin,...
                self.kwave,self.matrix * expansion.coefficients(:));
            
        end
        
        %-----------------------------------------
        % set the origin
        %-----------------------------------------

        function setOrigin(self,origin)
            
            self.origin = origin;
            
        end
        
    end % end methods
    
end

