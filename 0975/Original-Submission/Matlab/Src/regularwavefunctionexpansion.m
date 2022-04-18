% Regular wavefunction expansion
%
%  v = regularwavefunctionexpansion(u) copies a regular wavefunction 
%  expansion u into a new regular wavefunction expansion v.
%
%  v = regularwavefunctionexpansion(u,x) creates a new regular wavefunction
%  expansion v with origin x from the wavefunction expansion u using the
%  translation addition theorem to implement the change of origin.
%
%  v = regularwavefunctionexpansion(u,x,n) creates a new regular wavefunction
%  expansion v with order n and origin x from the wavefunction 
%  expansion u using the translation addition theorem to implement the change 
%  of origin.
%
%  v = regularwavefunctionexpansion(n,x,inc) where inc is of type 'incident' 
%  creates a wavefunction expansion v with order n, expansion origin x and 
%  coefficients taken from inc.
%
%  v = regularwavefunctionexpansion(n,x,k,cof) creates a wavefunction 
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
%  u = v + w returns regular wavefunction expansion obtained by adding
%  regular wavefunction expansions v and w.
%
%  u = v - w returns regular wavefunction expansion obtained by subtracting
%  regular wavefunction expansion w from regular wavefunction expansion v.
%
% Note: in the functions above vectors in the plane are represented by
% complex numbers.
%
% See also: wavefunctionexpansion, radiatingwavefunctionexpansion.
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


classdef regularwavefunctionexpansion < wavefunctionexpansion
    
    properties
        
    end
    
    methods
        
        %-----------------------------------------
        % constructor
        %-----------------------------------------
        
        function self = regularwavefunctionexpansion(varargin)
            
            % if we have one arguments then assume that varargin{1} is a
            % wavefunction expansion. Need to check that it is a regular
            % wavefunction expansion
            if nargin==1 && ~isa(varargin{1},'regularwavefunctionexpansion')
                error('Given wavefunction expansion must be a regularwavefunctionexpansion')
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
                
                if isa(varargin{1},'regularwavefunctionexpansion')
                    
                    if nargin == 2
                        % only change the origin
                        self.apply_addition_theorem(varargin{2},'SAME');
                    else
                        % change the origin and the order
                        self.apply_addition_theorem(varargin{2},'SAME',varargin{3});
                        self.order = varargin{3};
                    end
                    
                else
                    
                    if nargin == 2
                        % only change the origin
                        self.apply_addition_theorem(varargin{2},'DIFF');
                    
                    else
                        % change the origin and the order
                        self.apply_addition_theorem(varargin{2},'DIFF',varargin{3});
                        self.order = varargin{3};                        
                    end
                end
                
            end
            
        end
        
        %-----------------------------------------
        % add
        %-----------------------------------------
        
        function val = plus(self,other)
            
            % check types match
            if ~isa(other,'regularwavefunctionexpansion')
                error('Expansion types do not match')
            end
            
            % create a copy of self
            val = regularwavefunctionexpansion(self);
            
            % call parent add function
            plus@wavefunctionexpansion(val,other);
            
        end
        
        %-----------------------------------------
        % minus
        %-----------------------------------------
        
        function val = minus(self,other)
            
            % check types match
            if ~isa(other,'regularwavefunctionexpansion')
                error('Expansion types do not match')
            end
            
            % create a copy of self
            val = regularwavefunctionexpansion(self);
            
            % call parent add function
            minus@wavefunctionexpansion(val,other);
            
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
            
            val = sumcof(points,self.origin,self.kwave,self.coefficients,'J');
            
        end
        
    end
    
end

