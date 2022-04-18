% Incident field that is the difference of two other incident fields.
%
%   Note: this is intended to be used only by incident.minus method.
%
% See also: point_source, plane_wave, incident.
%
% Stuart C. Hawkins - 30 November 2018

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


classdef incidentminus < incident
    
    properties
        left
        right
    end
    
    methods
        
        %------------------------------------------------
        % constructor
        %------------------------------------------------
        
        function self = incidentminus(left,right)
        
            % check that input object is of class incident            
            if ~isa(left,'incident')
                error('left must be of class incident')
            end
            
            % check that input object is of class incident
            if ~isa(right,'incident')
                error('right must be of class incident')
            end
            
            % store given objects... we will use these when we need to
            % evaluate etc
            self.left = left;
            self.right = right;
            
        end
   
        %------------------------------------------------
        % get coefficients
        %------------------------------------------------

        function cof = get_coefficients(self,centre,nmax)
           
            % combine output of left and right object methods
            cof = self.left.get_coefficients(centre,nmax)...
                - self.right.get_coefficients(centre,nmax);
            
        end
    
        %------------------------------------------------
        % evaluate
        %------------------------------------------------

        function val = evaluate(self,varargin)
            
            % combine output of left and right object methods
            val = self.left.evaluate(varargin{:})...
                - self.right.evaluate(varargin{:});
        end
        
        %------------------------------------------------
        % evaluate gradient
        %------------------------------------------------

        function [dx,dy] = evaluateGradient(self,varargin)
           
            % use left and right object methods
            [left_dx,left_dy] = self.left.evaluateGradient(varargin{:});
            [right_dx,right_dy] = self.right.evaluateGradient(varargin{:});
            
            % combine output of left and right object methods
            dx = left_dx - right_dx;
            dy = left_dy - right_dy;
            
        end
                        
    end
    
end
    