% Plane wave object.
%
%  p = plane_wave(theta,k) returns a plane wave object p with 
%  wavenumber k and direction exp(1i*theta).
%
% Also:
%
%   f = p.evaluate(z) returns the values f of the plane wave at points z.
%
%   f = u.evaluate(z,mask) returns the values f of the plane wave at
%   points z for which mask==1 and NaN elsewhere.
%
%   [dx,dy] = u.evaluateGradient(z) returns dx and dy the partial 
%   derivatives of the plane wave in the x and y directions respectively
%   at the points z.
%
%   [dx,dy] = u.evaluateGradient(z,mask) returns dx and dy the partial 
%   derivatives of the plane wave in the x and y directions respectively
%   at the points z for which mask==1 and NaN elsewhere.
%
%   cof = p.get_coefficients(x,n) returns the vector cof of regular
%   wavefunction expansion coefficients of the plane wave field with 
%   wavefunction origin x and order n.
%
% Note: in the functions above vectors in the plane are represented by
% complex numbers.
%
% See also: point_source, incident.
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


classdef plane_wave < incident
   
    properties
        direction
        kwave
    end
    
    methods

        %-----------------------------------------------
        % constructor
        %-----------------------------------------------

        function self = plane_wave(angle,kwave)
            
            % set wavenumber 
            self.kwave = kwave;
            
            % set incident direction
            self.direction = exp(1i*angle);
            
        end
        
        %-----------------------------------------------
        % return vector of coefficients for scatterer at
        % centre
        %-----------------------------------------------

        function cof = get_coefficients(self,centre,nmax)
            
            % get the angle of the incident direction
            incang = angle(self.direction);
            
            n=-nmax:nmax;
            n=n(:);
            
            % set coefficients cf (3.13) in Ganesh, Hawkins, Hiptmair
            % IMA Journal of Numerical Analysis doi:10.1093/imanum/drr041
            cof = 1i.^abs(n).*exp(-1i*n*incang);

            % adjust the coefficients if the centre of the scatterer is not
            % the origin
            if centre~=0
                dp = real(conj(self.direction)*centre);
                phase = exp(1i*self.kwave*dp);
                cof = phase * cof;
            end
            
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

            % compute incident field
            v=exp(1i*self.kwave*real(conj(self.direction)*points));
            
            % insert values into the return array
            if nargin>2
                val(mask)=v;
            else
                val=v;
            end
            
        end
        
        %-----------------------------------------------
        % evaluate gradient
        %
        % This is useful for eg Neumann BCs. Note that we
        % cannot use complex numbers to represent vectors
        % in this case because the components of the vector
        % are, in general, complex.
        %-----------------------------------------------

        function [dx,dy] = evaluateGradient(self,points,mask)
            
            % initialize return values
            dx = zeros(size(points));
            dy = zeros(size(points));

            % apply mask if necessary
            if nargin>2
                points=points(mask);
            end

            % compute incident field
            scalarPart = 1i*self.kwave*exp(1i*self.kwave*real(conj(self.direction)*points));
            vx = real(self.direction) * scalarPart;
            vy = imag(self.direction) * scalarPart;
            
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