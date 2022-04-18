% Kite shaped obstacle
%
%  p = obstacleKite() creates an instance p of the obstacleKite class.
%
%  p.visualise() visualises the obstacle p.
%
%  [x,y,qx,qy,dqx,dqy,ddqx,ddqy,nqx,nqy,jac] = p.geom(t) computes points
%  (x,y) on the unit circle parametrised by t in [0,2*pi] and corresponding
%  points [qx,qy] on the kite shape. Also (dqx,dqy) are the derivatives of
%  qx and qy with respect to t. Also (ddqx,ddqy) are the second derivatives 
%  of qx and qy with respect to t. Also (nqx,nqy) is the unit outward
%  normal and jac is the Jacobian of the mapping from the unit circle to
%  the kite shape.
%
% See also: obstacle, obstaclePolar.
%
% See Colton and Kress, Inverse Acoustic and EM Scattering Theory 3rd Ed p79. 
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


classdef obstacleKite < obstacle
    
    properties        
    end
    
    methods
        
        %===============================================================
        % methods
        %===============================================================

        %-----------------------------------------
        % constructor
        %-----------------------------------------
        
        function self = obstacleKite()
           
            % call parent constructor
            self = self@obstacle();
                        
        end

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
        
        function [x,y,qx,qy,dqx,dqy,ddqx,ddqy,nqx,nqy,jac] = geom(self,t)
            
            % We use complex numbers to represent points in the plane. Then
            % z = x + 1i * y is associated with the point (x,y) in the
            % plane.

            % points on circle
            z = exp(1i*t);

            % get mapping and derivatives in complex form.. cf Colton and 
            % Kress, Inverse Acoustic and EM Scattering Theory 3rd Ed p79

            qz = cos(t) + 0.65*cos(2*t) - 0.65 + 1i * 1.5 * sin(t);            
            dqz = -sin(t) - 1.3*sin(2*t) + 1i * 1.5 * cos(t);            
            ddqz = -cos(t) - 2.6*cos(2*t) - 1i * 1.5 * sin(t);
            
            % now extract x and y parts
            x = real(z);
            qx = real(qz);
            dqx = real(dqz);
            ddqx = real(ddqz);
            y = imag(z);
            qy = imag(qz);
            dqy = imag(dqz);
            ddqy = imag(ddqz);
            
            % finally get jacobian
            jac = abs(dqz);
            nqx=dqy./jac;
            nqy=-dqx./jac;
            
        end
        
        
    end % end methods
    
end
        