% Radiating wavefunction expansion
%
%  v = radiatingwavefunctionexpansion(u) copies a radiating wavefunction
%  expansion u into a new radiating wavefunction expansion v.
%
%  v = radiatingwavefunctionexpansion(u,x) creates a new radiating wavefunction
%  expansion v with origin x from the radiating wavefunction expansion u using the
%  translation addition theorem to implement the change of origin.
%
%  v = radiatingwavefunctionexpansion(u,x,n) creates a new radiating wavefunction
%  expansion v with order n and origin x from the radiating wavefunction 
%  expansion u using the translation addition theorem to implement the change 
%  of origin.
%
%  v = radiatingwavefunctionexpansion(n,x,k,cof) creates a wavefunction
%  expansion v with order n, expansion origin x, wavenumber k and
%  coefficient vector cof.
%
% Also:
%
%  f = u.evaluate(z) returns the values f of the wavefunction expansion u
%  evaluated at the points z.
%
%  f = u.evaluate(z,mask) returns the values f of the wavefunction
%  expansion u evaluated at the points z for which mask==1 and NaN
%  elsewhere.
%
%  f = u.evaluateFarField(z) returns the values f of farfield of the
%  wavefunction expansion u evaluated at the points z on the unit circle,
%  ie abs(z) == 1.
%
%  u.changeorigin(x) changes the origin of the wavefunction expansion to x
%  using the translation addition theorem to update the expansion
%  coefficients.
%
%  u.rotatecoordinates(theta) converts the wavefunction expansion to a new
%  coordinate system obtained by rotation the old coordinate system by
%  angle theta in the positive direction. The expansion coefficients are
%  obtained accordingly.
%
%  c = u.getCoefficients() returns the vector of wavefunction expansion
%  coefficients.
%
%  u.visualize() visualizes the wavefunction expansion.
%
%  u.visualize(limits) visualizes the wavefunction expansion on the grid
%  with lower left corner limits(1) + 1i * limits(2) and upper right corner
%  limits(3) + 1i * limits(4).
%
%  u.visualize(limits,'-f') visualizes the wavefunction expansion with the
%  limit on the number of points overridden.
%
%  [f] = u.visualize(...) returns the value of the field f and does not
%  visualize it.
%
%  [f,z] = u.visualize(...) returns the value of the field f and the
%  corresponding points z and does not visualize it.
%
%  Warning: radiating wave functions may blow up in the plotting region.
%  Plotting manually using u.evaluate with a mask is recommended.
%
%  u.visualizeFarField() visualizes the far field of the wavefunction
%  expansion.
%
%  u.visualizeFarField(opts) visualizes the far field of the wavefunction
%  expansion with plotting options specified in opts.
%
%  u = v + w returns radiating wavefunction expansion obtained by adding
%  radiating wavefunction expansions v and w.
%
%  u = v - w returns radiating wavefunction expansion obtained by subtracting
%  radiating wavefunction expansion w from radiating wavefunction expansion v.
%
% Note: in the functions above vectors in the plane are represented by
% complex numbers.
%
% See also: wavefunctionexpansion, regularwavefunctionexpansion.
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


classdef radiatingwavefunctionexpansion < wavefunctionexpansion
    
    properties
        
        % none
        
    end
    
    methods
        
        %-----------------------------------------
        % constructor
        %-----------------------------------------
        
        function self = radiatingwavefunctionexpansion(varargin)
            
            % check that we aren't being given an incident field... the
            % parent constructor would work with this but the coefficients
            % would not be correct.
            if nargin>=3 && isa(varargin{3},'incident')
                error('Cannot get radiating field coefficients from incident field')
            end
            
            % if we have one arguments then assume that varargin{1} is a
            % wavefunction expansion. Need to check that it is a radiating
            % wavefunction expansion
            if nargin==1 && ~isa(varargin{1},'radiatingwavefunctionexpansion')
                error('Given wavefunction expansion must be a radiatingwavefunctionexpansion')
            end
            
            % call parent constructor
            self = self@wavefunctionexpansion(varargin{:});
            
            % if varargin{1} is a wavefunction expansion and we have more 
            % arguments then assume that varargin{2} is a new
            % origin... so we need to modify the coefficients using the
            % translation addition theorem
            % If a third argument is provided then the order is going to
            % be changed too.
            if isa(varargin{1},'wavefunctionexpansion') && nargin > 1
                
                if isa(varargin{1},'radiatingwavefunctionexpansion')
                    
                    if nargin == 2
                        % only change the origin
                        self.apply_addition_theorem(varargin{2},'SAME');
                    else
                        % change the origin and the order
                        self.apply_addition_theorem(varargin{2},'SAME',varargin{3});
                        self.order = varargin{3};
                    end
                    
                else
                    
                    error('Constructing radiating expansion from regular expansion is not supported.')
                    
                end
                
            end
                                   
        end
        
        %-----------------------------------------
        % evaluate far field
        %-----------------------------------------
        
        function val = evaluateFarField(self,points)
            
            % get far field values
            val = sumcof(points,self.origin,self.kwave,self.coefficients,'F');
            
        end
        
        %-----------------------------------------
        % visualize far field
        %-----------------------------------------

        % Just a wrapper for visualise.
        
        function varargout = visualizeFarField(varargin)
            
            [varargout{1:nargout}] = visualiseFarField(varargin{:});
            
        end
        
        %-----------------------------------------
        % visualise far field
        %-----------------------------------------
        
        function visualiseFarField(self,opts)
            
            % specify default options
            if nargin < 2
                opts = 'r-';
            end
            
            % set points
            n = 1000;
            t = 2*pi*(0:n-1)/n;
               
            % plot the far field
            plot(t,abs(self.evaluateFarField(exp(1i*t))),opts)
                
        end

        %-----------------------------------------
        % add
        %-----------------------------------------
        
        function val = plus(self,other)
            
            % check types match
            if ~isa(other,'radiatingwavefunctionexpansion')
                error('Expansion types do not match')
            end
            
            % create a copy of self
            val = radiatingwavefunctionexpansion(self);
            
            % call parent add function
            plus@wavefunctionexpansion(val,other);
            
        end
        
        %-----------------------------------------
        % minus
        %-----------------------------------------
        
        function val = minus(self,other)
            
            % check types match
            if ~isa(other,'radiatingwavefunctionexpansion')
                error('Expansion types do not match')
            end
            
            % create a copy of self
            val = radiatingwavefunctionexpansion(self);
            
            % call parent add function
            minus@wavefunctionexpansion(val,other);
            
        end
        
        %-----------------------------------------
        % uminus
        %-----------------------------------------
        
        function val = uminus(self)
            
            % copy
            val = self;
            
            % change the sign of the coefficients
            val.coefficients = -val.coefficients;
            
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
        
        function varargout = visualise(self,varargin)
                        
            if nargin < 2 || isempty(varargin{1})
            
                warning('Radiating wavefunctions may blow up in the visualisation region. Plotting manually using a mask is recommended.')
                
            end
            
            if nargin > 1 && ~isempty(varargin{1})
                
                % extract limits
                limits = varargin{1};
                
                % check if the origin is in the plot area
                if limits(1) < real(self.origin) && real(self.origin) < limits(2) ...
                        && limits(3) < imag(self.origin) && imag(self.origin) < limits(4)
                                       
                    warning('Radiating wavefunctions may blow up in the visualization region. Plotting manually using a mask is recommended.')
                
                end
            
            end
            
            % call parent method
            if nargout == 0
                visualise@wavefunctionexpansion(self,varargin{:});
            else
                varargout{:} = visualise@wavefunctionexpansion(self,varargin{:});                
            end
                
        end
        
    end % end methods
    
    methods(Access=protected)
        
        %-----------------------------------------
        % function to evaluate the expansion at given
        % points.
        % ** should be used only within methods
        % of the parent class. **
        %-----------------------------------------
        
        function val = internal_evaluate(self,points)
            
            val = sumcof(points,self.origin,self.kwave,self.coefficients,'H');
            
        end
        
    end
    
end

