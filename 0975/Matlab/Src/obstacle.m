% Obstacle
%
%  p = obstacle() creates an instance p of the obstacle class.
%
%  p.visualise() visualises the obstacle p.
%
%  Warning: this should be considered an ABSTRACT base class... it is not
%  intended that objects of this class be instantiated. Some methods of
%  this class must be overridden.
%
% See also: obstaclePolar, obstacleKite.
%
% Stuart C. Hawkins - 19 November 2014

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


classdef obstacle < handle
    
    properties        
    end
    
    methods
        
        %===============================================================
        % methods
        %===============================================================

        %-----------------------------------------
        % constructor
        %-----------------------------------------
        
        function self = obstacle()
           
            % nothing needed
            
        end
        
        %-----------------------------------------
        % visualise
        %-----------------------------------------
        
        function visualise(self,opts)
            
            % set default for opts
            if nargin<2
                opts = 'k-';
            end
            
            % evaluate the geometry parametrisation at lots of points and
            % plot
            t = 2*pi*(0:1000)/1000;
            [x,y,qx,qy] = self.geom(t);
            plot(qx,qy,opts)
            axis equal
            
        end
        
    end % end methods
    
    methods(Abstract=true)
        
        %===============================================================
        % these methods must be overridden in the child class
        %===============================================================

        %-----------------------------------------
        % geometry parametrisation
        %-----------------------------------------
        
        % (x,y) is a point on the circle
        % (qx,qy) is the corresponding point on the surface
        % (dqx,dqy) is the derivative wrt of (qx,qy)
        % (ddqx,ddqy) is the second derivative wrt of (qx,qy)
        % (nqx,nqy) is the normal to the surface at (qx,qy)
        % jac is the jacobian of the transformation from the circle to the
        % surface.
        
        [x,y,qx,qy,dqx,dqy,ddqx,ddqy,nqx,nqy,jac] = geom(self,t);
        
    end % end abstract methods
    
end
        