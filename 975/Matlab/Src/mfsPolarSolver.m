
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

classdef mfsPolarSolver < mfsSolver
    
    properties
        f
        df
    end
    
    methods
        
        %-----------------------------------------
        % constructor
        %-----------------------------------------
        
        function self = mfsPolarSolver(kwave,incidentField,f,df,n,tau,m)
            
            % set defaults for MFS parameters
            if nargin < 6
                tau = 5e-2;
            end
            
            if nargin < 7
                m = 2*n;
            end
            
            %  call parent constructor
            self = self@mfsSolver(kwave,incidentField,n,tau,m);
            
            % set radius function and its derivative
            self.f = f;
            self.df = df;
            
        end
                    
        %-----------------------------------------
        % visualize
        %-----------------------------------------
        
        % Just a wrapper for visualise.

        function varargout = visualize(varargin)

            [varargout{1:nargout}] = visualise(varargin{:});           
            
        end
        
        %-----------------------------------------
        % visualise
        %-----------------------------------------
        
        function visualise(self)
    
            theta = 2*pi*(0:1000)/1000;
            
            r = exp(1i*theta).*self.f(theta);
            
            plot(real(r),imag(r),'k-');
            
        end
        
    end % end methods
    
end