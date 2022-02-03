% Circular wavefunction
%
%  Warning: this should be considered an ABSTRACT base class... it is not
%  intended that objects of this class be instantiated. Several methods of
%  this class must be overridden.
%
% See also: regularwavefunction2d, radiatingwavefunction2d.
%
% Stuart C. Hawkins - 12 November 2014

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


classdef wavefunction2d < incident
    
    properties
        order
        kwave
        origin
    end
    
    methods

        %-------------------------------------------------
        % constructor
        %-------------------------------------------------
        
        function self = wavefunction2d(order,kwave,origin)

            % set default for origin
            if nargin < 3
                origin = 0;
            end
            
            % set properties
            self.order = order;
            self.kwave = kwave;
            self.origin = origin;
            
        end
        
        %-----------------------------------------------
        % function that determines radial behaviour
        %-----------------------------------------------

        % ** THIS FUNCTION MUST BE OVERLOADED IN THE CHILD CLASS **
        
        function val = radial_function(self,r)
            
            error('This function must be overloaded in child class')
            
        end
        
        %-----------------------------------------------
        % derivative of function that determines radial behaviour
        %-----------------------------------------------

        % ** THIS FUNCTION MUST BE OVERLOADED IN THE CHILD CLASS **

        function val = derivative_radial_function(self,r)
            
            error('This function must be overloaded in child class')
            
        end
        
        %-----------------------------------------------
        % evaluate method
        %-----------------------------------------------

        function val = evaluate(self,points,mask)
            
            % intialize return array
            val=zeros(size(points));

            % apply mask if necessary
            if nargin>2
                points=points(mask);
            end
            
            % subtract origin
            points = points - self.origin;

            % compute field
            v = exp(1i*self.order*angle(points)) ...
                .* self.radial_function(self.kwave * abs(points));
            
            % insert values into the return array
            if nargin>2
                val(mask)=v;
            else
                val=v;
            end
            
        end
        
        %-----------------------------------------------
        % evaluate method
        %-----------------------------------------------

        function [dx,dy] = evaluateGradient(self,points,mask)
            
            % intialize return array
            dx=zeros(size(points));
            dy=zeros(size(points));

            % apply mask if necessary
            if nargin>2
                points=points(mask);
            end
            
            % subtract origin
            points = points - self.origin;

            % wavefunction is given in polar coordinates so get necessary
            % quantities in polar coordinates first
            r = abs(points);
            f = self.radial_function(self.kwave * r);
            df = self.derivative_radial_function(self.kwave * r);            
            t = angle(points);
            dtdx = -imag(points)./r.^2;
            dtdy = real(points)./r.^2;
            
            % now compute gradient
            e = exp(1i*self.order*t);
            vx = 1i*self.order*e.*f.*dtdx + self.kwave*e.*df.*real(points)./r;
            vy = 1i*self.order*e.*f.*dtdy + self.kwave*e.*df.*imag(points)./r;
            
            % insert values into the return array
            if nargin>2
                dx(mask) = vx;
                dy(mask) = vy;
            else
                dx = vx;
                dy = vy;
            end
            
        end
        
    end % end methods
    
end
