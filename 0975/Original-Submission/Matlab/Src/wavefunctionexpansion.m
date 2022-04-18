% Wavefunction expansion
%
%  v = wavefunctionexpansion(u) copies a wavefunction expansion u into
%  a new wavefunction expansion v.
%
%  v = wavefunctionexpansion(n,x,inc) where inc is of type 'incident' 
%  creates a wavefunction expansion v with order n, expansion origin x and 
%  coefficients taken from inc.
%
%  v = wavefunctionexpansion(n,x,k,cof) creates a wavefunction expansion v 
%  with order n, expansion origin x, wavenumber k and coefficient vector 
%  cof.
%
% Also:
%
%  c = v.getCoefficients() returns the vector of wavefunction expansion
%  coefficients.
%
%  Warning: this should be considered an ABSTRACT base class... it is not
%  intended that objects of this class be instantiated. Several methods of
%  this class must be overridden.
%
% See also: regularwavefunctionexpansion, radiatingwavefunctionexpansion.
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


classdef wavefunctionexpansion < handle
    
    properties
        order
        origin
        kwave
        coefficients
    end
    
    methods
        
        %-----------------------------------------
        % constructor
        %-----------------------------------------
        
        function self = wavefunctionexpansion(varargin)
            
            if (nargin==1 || nargin==2 || nargin==3) && isa(varargin{1},'wavefunctionexpansion')
                                
                self.order = varargin{1}.order;
                self.origin = varargin{1}.origin;
                self.kwave = varargin{1}.kwave;
                self.coefficients = varargin{1}.coefficients(:);
                                
            elseif nargin==3
                
                self.order = varargin{1};
                self.origin = varargin{2};
                
                if ~isa(varargin{3},'incident')
                    
                    error('Third argument must be an incident wave')
                    
                end
                
                self.kwave = varargin{3}.kwave;
                
                self.coefficients = varargin{3}.get_coefficients(self.origin,...
                    self.order);
                
            elseif nargin==4
                
                % check that the coefficient vector is of the correct size
                if length(varargin{4}) ~= 2*varargin{1}+1
                    
                    error('The order n does not match the size of the coefficient vector cof.')
                    
                end
                
                self.order = varargin{1};
                self.origin = varargin{2};
                
                self.kwave = varargin{3};
                self.coefficients = varargin{4}(:);
                            
            else
                
                error('Incorrect parameters to constructor')
                
            end
            
            
        end
        
        %-----------------------------------------
        % evaluate near field
        %-----------------------------------------
        
        function val = evaluate(self,points,mask)
            
            % initialize array
            val = zeros(size(points));
            
            % apply mask if necessary
            if nargin>2
                points=points(mask);
            end
            
            % get far field values
            temp = self.internal_evaluate(points);
            
            % insert into return array
            if nargin>2
                val(mask) = temp;
                val(~mask) = NaN;
            else
                val = temp;
            end
            
        end
        
        %-----------------------------------------
        % change origin
        %-----------------------------------------
        
        function changeorigin(self,new_origin)
            
            % call addition theorem function
            self.apply_addition_theorem(new_origin,'SAME');
            
        end
        
        %-----------------------------------------
        % rotate coordinate system
        %-----------------------------------------

        function varargout = rotatecoordinates(self,angle)
            
            % compute diagonal of rotation matrix
            alpha = exp(1i*(-self.order:self.order)*angle);
            
            % apply rotation matrix (which is diagonal)
            self.coefficients = alpha(:).* self.coefficients;
            
            % return self if necessary
            if nargout > 0
                varargout{1} = self;
            end
            
        end

        %-----------------------------------------
        % add
        %-----------------------------------------
        
        function plus(self,other)
            
            % check that the parameters match
            if self.order ~= other.order
                error('Orders do not match')
            end
            if self.origin ~= other.origin
                error('Origins do not match')
            end
            if self.kwave ~= other.kwave
                error('Wavenumbers do not match')
            end
            
            % add the coefficients
            self.coefficients = self.coefficients + other.coefficients;
            
        end
        
        %-----------------------------------------
        % minus
        %-----------------------------------------
        
        function minus(self,other)
            
            % check that the parameters match
            if self.order ~= other.order
                error('Orders do not match')
            end
            if self.origin ~= other.origin
                error('Origins do not match')
            end
            if self.kwave ~= other.kwave
                error('Wavenumbers do not match')
            end
            
            % add the coefficients
            self.coefficients = self.coefficients - other.coefficients;
            
        end
        
        %-----------------------------------------
        % get coefficients
        %-----------------------------------------

        function val = getCoefficients(self)
            
            val = self.coefficients;
            
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

        function varargout = visualise(self,limits,opts)                       

            % - - - - - - - - - - - - - - -
            % set some parameters
            % - - - - - - - - - - - - - - -

            % points per wavelength
            ppw = 15;
            
            % maximum number of points (without override)
            maxpts = 300^2;
               
            % default radius in terms of wavelength
            dr = 4;
            
            % - - - - - - - - - - - - - - -
            % set some parameters
            % - - - - - - - - - - - - - - -

            % determine wavelength
            lambda = 2*pi / self.kwave;
            
            % - - - - - - - - - - - - - - -
            % set limits
            % - - - - - - - - - - - - - - -
            
            % set default for axis limits if none provided
            if nargin < 2 || isempty(limits)
                limits = [real(self.origin) - dr * lambda, real(self.origin) + dr * lambda, ...
                    imag(self.origin) - dr * lambda, imag(self.origin) + dr * lambda];
            end
                        
            % determine number of points required
            nx = ceil( ppw * sqrt(2) * (limits(2)-limits(1)) / lambda);
            ny = ceil( ppw * sqrt(2) * (limits(4)-limits(3)) / lambda);
                        
            % apply maximum to number of points
            if nx*ny > maxpts && (nargin<3 || ~strcmp(opts,'-f'))
                
                % display warning
                warning('Requested number of points %d is being restricted to %d. Use opts = ''-f'' to force. The CPU time may be long.',nx*ny,maxpts);
                
                % limit number of points
                tmp = sqrt(maxpts/(nx*ny));
                nx = ceil( tmp * nx );
                ny = ceil( tmp * ny );
                
            end
                        
            % - - - - - - - - - - - - - - -
            % setup grid and plot
            % - - - - - - - - - - - - - - -

            % setup a grid
            s=linspace(limits(1),limits(2),nx);
            t=linspace(limits(3),limits(4),ny);
            [x,y]=meshgrid(s,t);
            z = x+y*1i;            

            % get the field
            u = self.evaluate(z);
            
            % plot
            if nargout == 0
                surf(x,y,real(u));
                view([0 90]);
                shading interp;
                colorbar
            end
            
            % - - - - - - - - - - - - - - -
            % setup return values
            % - - - - - - - - - - - - - - -

            if nargout > 0
                
                % return the computed field
                varargout{1} = u;
                
            end
            
            if nargout > 1
                
                % return the points used
                varargout{2} = z;
                
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
        
        % ** must be overloaded in child class **
        
        function val = internal_evaluate(self,points)
            
            error('Not implemented yet')
            
        end
                
        %-----------------------------------------
        % function to apply the addition theorem
        %-----------------------------------------
        
        % See (36)--(37) in M. Ganesh and S. C. Hawkins, JQSRT 123
        % (2013) 41--52.
        % This part can be computed quickly once the transfer matrix
        % self.transfer_matrix has been computed and stored.
        
        function apply_addition_theorem(self,new_origin,type,new_order)
                        
            % set default for new_order
            if nargin<4
                new_order = self.order;
            end
            
            % get change of origin vector
            vec = new_origin - self.origin;
            
            % get polar coordinates
            r = abs(vec);
            theta = angle(vec);
            
            % setup arrays of indices
            [nn,mm]=meshgrid(-self.order:self.order,-new_order:new_order);
            pp=nn-mm;
            
            % get Bessel function
            if strcmp(type,'SAME')
                bess = besselj(0:2*self.order,self.kwave*r);
            else
                bess = besselh(0:2*self.order,self.kwave*r);
            end
            B=bess(abs(pp)+1);
            
            % vectorized way to compute the transfer matrix
            matrix = B.*exp(1i*pp*theta).*wavefunctionexpansion.get_transfer_matrix(new_order,self.order);
            
            % update the coefficients
            self.coefficients = matrix * self.coefficients;
            
            % update the origin
            self.origin = new_origin;
            
        end
        
    end % end methods
    
    methods(Static)
        
        %-----------------------------------------
        % function to compute the transfer matrix
        %-----------------------------------------
        
        % See (39) in M. Ganesh and S. C. Hawkins, JQSRT 123
        % (2013) 41--52. This part is independent of the translation
        % and hence can be precomputed and stored.
        %
        % Here we use a static method and persistent variables to store
        % the transfer matrix once and re-use it for all wavefunction
        % expansion objects. If a transfer matrix of higher order is
        % required then the higher order transfer matrix is stored. If a 
        % transfer matrix of lower order is subsequently required then the
        % appropriate subset is used. 
        % Note: this may lead to inconsistent accuracy but this is worth it 
        % for speed and reduced storage.
        
        function varargout = get_transfer_matrix(orderout,orderin)
        
            % declare persistent variables... these will behave a bit like 
            % static variables in other OO languages
            persistent nmax
            persistent transfer_matrix            

            % get order as larger of orderin and orderout
            order = max(orderin,orderout);
            
            % we use nmax being empty as a sign that the transfer matrix
            % has not been created yet.
            if isempty(nmax)
                
                if nargin < 1
                    error('order must be provided');
                else 
                    
                    % initialize nmax and transfer_matrix
                    nmax = 0;
                    transfer_matrix = [];
                
                end
                
            end
            
            % we assume we have the nmax transfer matrix. If the desired
            % order is bigger than nmax then we will need to recompute the
            % transfer matrix.
            if order > nmax
            
                % for brevity denote order by N
                N = order;
                
                % setup quadrature points
                tp = 2*pi*(0:4*N+3)/(4*N+4);
                tw = 2*pi/(4*N+4);
                
                % initialize matrix
                transfer_matrix = zeros(2*N+1,2*N+1);
                                
                for n=-N:N
                    for m=-N:N
                        
                        p = n-m;
                        
                        % See (39) in M. Ganesh and S. C. Hawkins, JQSRT 123
                        % (2013) 41--52
                        transfer_matrix(m+N+1,n+N+1) = ...
                            sqrt(2*pi)*(-1i).^(abs(n)-abs(m)-abs(p)) ...
                            *sum(tw.*exp(-1i*m*tp).*exp(1i*n*tp) ...
                            .*exp(-1i*p*tp))/sqrt(2*pi)^3;
                        
                    end
                end
            
                % store the size of the computed transfer matrix
                nmax = order;
                
            end

            % return the transfer matrix (or appropriate portion of it)
            if nargout > 0
                varargout{1} = transfer_matrix((-orderout:orderout)+nmax+1,(-orderin:orderin)+nmax+1);
            end
            
        end
        
    end
    
end