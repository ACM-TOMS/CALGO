% Incident field
%
%  Warning: this should be considered an ABSTRACT base class... it is not
%  intended that objects of this class be instantiated. Several methods of
%  this class must be overridden.
%
% See also: plane_wave, point_source, regularwavefunction2d,
% radiatingwavefunction2d.
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


classdef incident < handle
    
    properties
    end
    
    methods
        
        % This is an abstract base class for the plane wave and 
        % point source type incident fields. The methods
        % just do the minimum... the idea is that they
        % get replaces by specific versions for the other
        % incident fields.
        
        %-----------------------------------------------
        % constructor
        %-----------------------------------------------

        function cof = get_coefficients(self,centre,nmax)
            
            cof=zeros(2*nmax+1,1);
            
        end
        
        
        %-----------------------------------------------
        % evaluate method
        %-----------------------------------------------

        function val = evaluate(self,points,mask)
            
            val=zeros(size(points));
            
        end

        %-----------------------------------------------
        % evaluate method
        %-----------------------------------------------

        function [dx,dy] = evaluateGradient(self,points,mask)
            
            error('Not implemented yet')
            
            dx = zeros(size(points));
            dy = zeros(size(points));
            
        end
        
    end
    
end