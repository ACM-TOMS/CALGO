% Obstacle specified by polar coordinates.
%
%  p = obstacle(r,dr,ddr) creates an instance p of the obstacle given in
%  polar coordinates where r,dr,ddr are anonymous functions specifying the 
%  radius and its first and second derivatives with respect to the angle.
%
%  p.visualise() visualises the obstacle p.
%
%  [x,y,qx,qy,dqx,dqy,ddqx,ddqy,nqx,nqy,jac] = p.geom(t) computes points
%  (x,y) on the unit circle parametrised by t in [0,2*pi] and corresponding
%  points [qx,qy] on the obstacle. Also (dqx,dqy) are the derivatives of
%  qx and qy with respect to t. Also (ddqx,ddqy) are the second derivatives 
%  of qx and qy with respect to t. Also (nqx,nqy) is the unit outward
%  normal and jac is the Jacobian of the mapping from the unit circle to
%  the obstacle shape.
%
% Example: trefoil shape
%
%  r = @(t) 1 + 0.3 * cos(3*t);
%  dr = @(t) -0.9*sin(3*t);            
%  ddr = @(t) -2.7*cos(3*t);            
%  p = obstacle(r,dr,ddr);
%  p.visualise()
%  
% See also: obstacle, obstacleKite.
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


classdef obstaclePolar < obstacle
    
    properties        
        r
        dr
        ddr
    end
    
    methods
        
        %===============================================================
        % methods
        %===============================================================

        %-----------------------------------------
        % constructor
        %-----------------------------------------
        
        function self = obstaclePolar(r,dr,ddr)
           
            % call parent constructor
            self = self@obstacle();
            
            % set radius and its derivatives
            self.r = r;
            self.dr = dr;
            self.ddr = ddr;
            
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
            
            % radius etc
            r = self.r(t);
            dr = self.dr(t);
            ddr = self.ddr(t);
            
            % qmap etc
            qz = r .* z;
            dqz = (dr + 1i * r) .* z;
            ddqz = ( ddr + 2i * dr - r ) .* z;
            
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
        