
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

classdef mfsPolarSoftSolver < mfsPolarSolver
    
    properties
    end
    
    methods
        
        %-----------------------------------------
        % constructor
        %-----------------------------------------
        
        function self = mfsPolarSoftSolver(kwave,incidentField,f,df,n,tau,m)
            
            % set defaults for MFS parameters
            if nargin < 6
                tau = 5e-2;
            end
            
            if nargin < 7
                m = 2*n;
            end
            
            %  call parent constructor
            self = self@mfsPolarSolver(kwave,incidentField,f,df,n,tau,m);            
        
        end
            
        %-----------------------------------------
        % setup
        %-----------------------------------------
        
        function setup(self)

            % setup structure with the parameters for the MFS method
            opts= struct('eta',self.kwave,'fast',2,'multiplier',2.1,'tau',self.tau);
            
            % create boundary segment from the parametrisation
            boundary = segment.radialfunc(self.m, {self.f,self.df});
            
            % set a homogeneous Neumann boundary condition on the boundary
            % ie sound-soft BC
            boundary.setbc(1,'D', []);
            
            % setup a domain outside the boundary
            d = domain([], [], boundary, -1);
            
            % setup a MFS basis on the boundary
            d.addmfsbasis(boundary,self.n,opts);
            
            % initialise the scattering problem
            self.scatteringObject = scattering(d, []);

        end
        
    end % end methods
    
end