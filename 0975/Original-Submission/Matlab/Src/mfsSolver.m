
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

classdef mfsSolver < solver
    
    properties
        scatteringObject
        tau
        m
        n
        coeffs
    end
    
    methods
        
        %-----------------------------------------
        % constructor
        %-----------------------------------------
        
        function self = mfsSolver(kwave,incidentField,n,tau,m)
            
            %  call parent constructor
            self = self@solver(kwave,incidentField)
            
            % set defaults for MFS parameters
            if nargin < 4
                tau = 5e-2;
            end
            
            if nargin < 5
                m = 2*n;
            end
            
            % set key parameters
            self.m = m;
            self.n = n;
            self.tau = tau;

            % set derived parameters to empty... these will be created when
            % setup is called
            self.scatteringObject = [];
            self.coeffs = [];
            
        end
        
        %-----------------------------------------
        % solve
        %-----------------------------------------
        
        function solve(self)
            
            if isempty(self.scatteringObject)
                
                error('Must call setup() first')
                
            end
            
            % set wavenumber
            self.scatteringObject.setoverallwavenumber(self.kwave);

            %- - - - - - - - - - - - - - - - - - -
            % setup right hand side and solve
            %- - - - - - - - - - - - - - - - - - -
            
            % here we manually do some of the functionality of bvp.solvecoeffs
            % so that it can handle multiple right hand sides... we setup
            % the RHS to be a matrix and then the standard bvp.solvecoeffs
            % code will solve for all RHS simultaneously.
            
            % set the quadrature weights
            self.scatteringObject.fillquadwei();
            
            % loop through incident directions
            for k=1:length(self.incidentField)
                
                % set incident wave direction
                if isa(self.incidentField{k},'plane_wave')
                    
                    % set the incident wave in MPSPACK
                    self.scatteringObject.setincidentwave(angle(self.incidentField{k}.direction));
                                        
                elseif isa(self.incidentField{k},'regularwavefunction2d')

                    % define functions for the incident wave and its
                    % derivatives
                    ui = @(x) self.incidentField{k}.evaluate(x);
                    uix = @(x) self.incidentField{k}.evaluateGradient(x);
                    f = @(x) self.incidentField{k}.evaluateGradient(x);
                    uiy = @(x) getSecondOutput(f,x);
                    
                    % set the incident wave in MPSPACK
                    self.scatteringObject.setincidentwave(ui,uix,uiy);
                    
                else
                    
                    error('incident field must be a plane wave or a regularwavefunction2d')
                    
                end
            
                self.scatteringObject.rhs(:,k) = self.scatteringObject.fillrighthandside();
                
            end
                            
            % solve scattering problem
            self.scatteringObject.solvecoeffs;
                
            for k=1:length(self.incidentField)

                % store coefficients
                self.coeffs{k} = self.scatteringObject.co(:,k);
                
            end
            
        end
        
        %-----------------------------------------
        % get far field
        %-----------------------------------------
        
        function val = getFarField(self,points,index)
            
            if nargin<3
                index = 1;
            end
            
            if isempty(self.coeffs)
                
                error('Must run solve() first')
                
            end
            
            for k=1:length(index)
            
                % set coefficients
                self.scatteringObject.co = self.coeffs{index(k)};
                
                % get the far field
                opts = [];
                val(:,k) = self.scatteringObject.gridfarfield(opts,points);
            
            end
                
        end

        %-----------------------------------------
        % get field
        %-----------------------------------------
        
        function val = getField(self,points,index)

            if nargin<3
                index = 1;
            end
            
            if isempty(self.coeffs)
                
                error('Must run solve() first')
                
            end
            
            % get the size of points
            [n,m] = size(points);
            
            ii=abs(points)>1/sqrt(3);

            % reshape points into a vector
            points = reshape(points(find(ii)),[],1);
            
            p=pointset(points);
            
            for k=1:length(index)
            
                % set coefficients
                self.scatteringObject.co = self.coeffs{index(k)};
                
                % get the far field
                opts = [];
                val(ii,k) = self.scatteringObject.pointsolution(p);
            
            end

            % reshape val into the original shape of points
            val = reshape(val,n,m,[]);
            
        end
    
    end % end methods
    
end