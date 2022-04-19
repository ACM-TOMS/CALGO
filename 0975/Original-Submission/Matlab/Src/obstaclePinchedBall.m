% Pinched ball shaped obstacle
%
%  p = obstaclePinchedBall() creates an instance p of the 
%  obstaclePinchedBall class.
%
%  p.visualise() visualises the obstacle p.
%
%  [x,y,qx,qy,dqx,dqy,ddqx,ddqy,nqx,nqy,jac] = p.geom(t) computes points
%  (x,y) on the unit circle parametrised by t in [0,2*pi] and corresponding
%  points [qx,qy] on the pinched ball. Also (dqx,dqy) are the derivatives of
%  qx and qy with respect to t. Also (ddqx,ddqy) are the second derivatives 
%  of qx and qy with respect to t. Also (nqx,nqy) is the unit outward
%  normal and jac is the Jacobian of the mapping from the unit circle to
%  the pinched ball.
%
% See also: obstacle, obstaclePolar.
%
% See Colton and Kress, Inverse Acoustic and EM Scattering Theory 3rd Ed p94. 
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


classdef obstaclePinchedBall < obstaclePolar
    
    properties
    end
    
    methods
        
        %-----------------------------------------
        % constructor
        %-----------------------------------------
        
        function self = obstaclePinchedBall()
           
            % set anonymous functions cf Colton and Kress, Inverse Acoustic 
            % and EM Scattering Theory 3rd Ed p94 with theta = pi/2
            r = @(t) sqrt(1.44 - 0.5*cos(-4*t));
            dr = @(t) -sin(-4*t)./sqrt(1.44 - 0.5*cos(-4*t));
            ddr = @(t) -sin(-4*t).^2./(1.44 - 0.5*cos(-4*t)).^1.5 ...
                + 4*cos(-4*t)./sqrt(1.44 - 0.5*cos(-4*t));
            
            % call parent constructor
            self = self@obstaclePolar(r,dr,ddr);
                        
        end
        
    end
    
end