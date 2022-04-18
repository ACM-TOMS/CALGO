
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

classdef solver_name < solver

    % ** Note: the class name above should match the filename and the name
    % of the constructor below **
    
    properties
        % ** put properties required for your solver here **
    end
    
    methods
        

        %-----------------------------------------
        % constructor
        %-----------------------------------------
        
        function self = solver_name(kwave,incidentField)

            % ** Note: the solver_name above should be replaced with the
            % name of your class. You may add other parameters as required.
            % **
            
            % call parent constructor
            self = self@solver(kwave,incidentField)
            
            % ** code to initialize your solver class goes here **
            
        end
                 
        %===============================================================
        % these methods must be provided 
        %===============================================================

        %-----------------------------------------
        % setup
        %-----------------------------------------
        
        % This methods sets up your solver, eg assembles discretisation
        % matrices etc
        
        function setup(self)                      
            
            % ** your code goes here **
            
        end
        
        %-----------------------------------------
        % solve
        %-----------------------------------------
        
        % This method solves the scattering problem for every right hand
        % side specified in the self.incidentField cell array.
        
        function solve(self)

            % ** recommend you check here that the setup method has been
            % run **
            
            % ** your code goes here **            
            
        end
        
        %-----------------------------------------
        % get far field
        %-----------------------------------------
        
        % This method should compute the far field for the incident fields
        % self.incidentField{k} for each k in the array index. The return
        % value val should be an array with the column val(:,j) containing
        % the far field for self.incidentField{index(j)}.
        
        function val = getFarField(self,points,index)
            
            % ** your code goes here **            

        end
        
        %===============================================================
        % you may provide other methods required to implement your solver
        % or help the user
        %===============================================================

    end % end methods
    
end