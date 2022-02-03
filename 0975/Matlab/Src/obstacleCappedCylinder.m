% Capped cylinder shaped obstacle
%
%  p = obstacleCappedCylinder(h,w) creates an instance p of the 
%  obstacleCappedCylinder class. The cylinder body has height h and width w.
%  The total height of the cylinder is h+w, including the two semi-circular
%  caps with radius w/2.
%
%  p.visualise() visualises the obstacle p.
%
%  [x,y,qx,qy,dqx,dqy,ddqx,ddqy,nqx,nqy,jac] = p.geom(t) computes points
%  (x,y) on the unit circle parametrised by t in [0,2*pi] and corresponding
%  points [qx,qy] on the capped cylinder. Also (dqx,dqy) are the derivatives of
%  qx and qy with respect to t. Also (ddqx,ddqy) are the second derivatives 
%  of qx and qy with respect to t. Also (nqx,nqy) is the unit outward
%  normal and jac is the Jacobian of the mapping from the unit circle to
%  the capped cylinder.
%
% See also: obstacle, obstaclePolar, obstacleKite.
%
% Stuart C. Hawkins - 16 November 2015

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

classdef obstacleCappedCylinder < obstacle
    
    properties        
        height
        width
    end
    
    methods
        
        %===============================================================
        % methods
        %===============================================================

        %-----------------------------------------
        % constructor
        %-----------------------------------------
        
        function self = obstacleCappedCylinder(height,width)
           
            % call parent constructor
            self = self@obstacle();
                        
            % set properties
            self.height = height;
            self.width = width;
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

            % get mapping and derivatives in complex form.. 
            
            % compute related quantities
            r = self.width / 2;
            d = self.height / 2;
            hp = self.height + pi * r;
            c = hp/pi;
            
            % map t to arc length s
            s = c * t;
            
            % compute segment boundaries
            s1 = self.height;
            s2 = self.height + r*pi;
            s3 = 2*self.height + r*pi;         
            
            % initialise arrays
            qz = zeros(size(t));
            dqz = zeros(size(t));
            ddqz = zeros(size(t));
            
            % first segment - right edge
            ii = (s <= s1);
            qz(ii) = r + 1i*(s(ii)-d);
            dqz(ii) = c*1i;
            ddqz(ii) = zeros(size(s(ii)));
            
            % second segment - top circle
            ii = (s1 < s) & (s < s2);
            qz(ii) = d*1i + r * exp(1i*(s(ii)-s1)/r);
            dqz(ii) = c*1i * exp(1i*(s(ii)-s1)/r);
            ddqz(ii) = -c^2*(1/r) * exp(1i*(s(ii)-s1)/r);
            
            % third segment - left edge
            ii = (s2 <= s) &  ( s <= s3);
            qz(ii) = -r + 1i*(d-(s(ii)-s2));
            dqz(ii) = -c*1i;
            ddqz(ii) = zeros(size(s(ii)));
            
            % last segment - bottom circle
            ii = (s3 < s);
            qz(ii) = - d*1i + r * exp(1i*(pi+(s(ii)-s3)/r));
            dqz(ii) = c*1i * exp(1i*(pi+(s(ii)-s3)/r));
            ddqz(ii) = -c^2*(1/r) * exp(1i*(pi+(s(ii)-s3)/r));

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
        