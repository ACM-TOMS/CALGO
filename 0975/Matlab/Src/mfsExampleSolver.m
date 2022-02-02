
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

classdef mfsExampleSolver < solver
    
    properties
        f
        df
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
        
        function self = mfsExampleSolver(kwave,incidentField,f,df,n,tau,m)
            
            % set defaults for MFS parameters
            if nargin < 6
                tau = 5e-2;
            end
            
            if nargin < 7
                m = 2*n;
            end
            
            %  call parent constructor
            self = self@solver(kwave,incidentField);
            
            % set radius function and its derivative
            self.f = f;
            self.df = df;
            
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

            % set wavenumber
            self.scatteringObject.setoverallwavenumber(self.kwave);

        end
        
        %-----------------------------------------
        % solve
        %-----------------------------------------
        
        function solve(self)
            
            if isempty(self.scatteringObject)
                
                error('Must call setup() first')
                
            end
            
            %- - - - - - - - - - - - - - - - - - -
            % setup right hand side and solve
            %- - - - - - - - - - - - - - - - - - -
            
            % loop through incident directions
            for k=1:length(self.incidentField)
                
                % define functions for the incident wave and its
                % derivatives
                ui = @(x) self.incidentField{k}.evaluate(x);
                uix = @(x) self.incidentField{k}.evaluateGradient(x);
                f = @(x) self.incidentField{k}.evaluateGradient(x);
                uiy = @(x) getSecondOutput(f,x);
                
                % set the incident wave in MPSPACK
                self.scatteringObject.setincidentwave(ui,uix,uiy);
                
                % clear the RHS
                self.scatteringObject.rhs = [];
                
                % solve scattering problem
                self.scatteringObject.solvecoeffs;
                
                % store coefficients
                self.coeffs{k} = self.scatteringObject.co;
                
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
        % visualise
        %-----------------------------------------
        
        function visualise(self)
    
            theta = 2*pi*(0:1000)/1000;
            
            r = exp(1i*theta).*self.f(theta);
            
            plot(real(r),imag(r),'k-');
            
        end
        
    end % end methods
    
end