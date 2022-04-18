% Translation addition theorem matrix.
%
%   obj = translation(n) sets up an instance of class translation which
%   implements the translation addition theorem for wavefunction series 
%   expansions with order up to n.
%
%   M = obj.translationMatrixOutToIn(z,waveNumber,n) returns n x n
%   translation addition theorem matrix M. 
%
%   If b is a vector of radiating wavefunction expansion coefficients then 
%   a = M * b is a vector of regular wavefunction expansion coefficients 
%   with origin translated by z. Here waveNumber is the wavenumber of the 
%   series expansions.
%
%   M = obj.translationMatrixOutToOut(z,waveNumber,n) returns n x n
%   translation addition theorem matrix M. 
%
%   If b is a vector of radiating wavefunction expansion coefficients then 
%   a = M * b is a vector of radiating wavefunction expansion coefficients 
%   with origin translated by z. Here waveNumber is the wavenumber of the 
%   series expansions.
%
%   If b is a vector of regular wavefunction expansion coefficients then 
%   a = M * b is a vector of regular wavefunction expansion coefficients 
%   with origin translated by z. Here waveNumber is the wavenumber of the 
%   series expansions.
%
%   The translation addition theorem coefficients are given by (37)--(39)
%   of [1] after accounting for the different normalisation of the
%   wavefunctions, or derived using identical arguments used in [2] for
%   spherical wavefunctions.
% 
%   References:
%
%   [1] M. Ganesh and S. C. Hawkins. J. Quant. Spectrosc. Radiat. Transfer
%   123 (2013), 41--52.
%
%   [2] T. J. Dufva et al. Progress in Electromagnetics Research B 4
%   (2008), 79--99.
%
% Note: in the above vectors in the plane are represented by
% complex numbers.
%
% Stuart C. Hawkins - 6 December 2018

% Copyright 2014-2019 Stuart C. Hawkins
% 	
% This file is part of MIESOLVER.
% 
% MIESOLVER is free software: you can redistribute it and/or modify	
% it under the terms of the GNU General Public License as published by	
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% MIESOLVER is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with MIESOLVER.  If not, see <http://www.gnu.org/licenses/>.


classdef translation < handle
    
    properties
        nmax
        transferMatrix
    end
    
    methods
        
        %------------------------------------------------------------------
        % contructor: this also sets up the transfer matrix, which is used
        % to compute all transfer matrices, regardless of wavenumber of
        % translation vector
        %------------------------------------------------------------------
        
        function self = translation(nmax)
            
            %----------------------------
            % preliminaries
            %----------------------------
            
            % set nmax
            self.nmax = nmax;
            
            %----------------------------
            % create transfer matrix
            %----------------------------

            % set up quadrature points and weights
            tp = 2*pi*(0:4*self.nmax+3)/(4*self.nmax+4);
            tw = 2*pi/(4*self.nmax+4);
            
            % initialise matrix
            self.transferMatrix = zeros(2*self.nmax+1,2*self.nmax+1);
            
            % compute matrix entries
            for n=-self.nmax:self.nmax
                for m=-self.nmax:self.nmax
                    
                    p = n-m;
                    
                    self.transferMatrix(m+self.nmax+1,n+self.nmax+1) = ...
                        sqrt(2*pi)*(-1i).^(abs(n)-abs(m)-abs(p)) ...
                        *sum(tw.*exp(-1i*m*tp).*exp(1i*n*tp).*exp(-1i*p*tp))...
                        /sqrt(2*pi)^3;
                    
                end
            end
            
        end
           
        %------------------------------------------------------------------
        % compute translation matrix to write a radiating field as an
        % incident field
        %------------------------------------------------------------------

        function matrix = translationMatrixOutToIn(self,translationVector,waveNumber,n)
            
            % Note: translation vector is given as a complex number with
            % real and imaginary parts corresponding to the x and y
            % coordinates respectively
                       
            % process translation vector
            r = abs(translationVector);
            theta=angle(translationVector);

            % in the addition theorem we only sum over some indices so
            % determine those indices here; ii holds the indices that we
            % sum over
            [nn,mm]=meshgrid(-self.nmax:self.nmax,-self.nmax:self.nmax);
            pp=nn-mm;

            % initialise matrix
            matrix = zeros(2*self.nmax+1,2*self.nmax+1);

            % get the necessary Hankel function values
            bess = besselh(0:2*self.nmax,waveNumber*r);
            B=bess(abs(pp)+1);

            % vectorised way to compute the transfer matrix
            matrix = B.*exp(1i*pp*theta).*self.transferMatrix;
            
            % process parameters
            if nargin > 3
                
                if max(n) > self.nmax
                    error('n must not be greater than self.nmax')
                end
                
                % setup indices into matrix
                indexi = (-n(1):n(1)) + self.nmax + 1;
                indexj = (-n(2):n(2)) + self.nmax + 1;
                
                % return only the relevant part of the matrix
                matrix = matrix(indexi,indexj);
                
            end
            
        end
        
        %------------------------------------------------------------------
        % compute translation matrix to write a radiating field as an
        % radiating field with a different origin. Same matrix writes an
        % incident field as an incident field with a different origin
        %------------------------------------------------------------------
        
        function matrix = translationMatrixOutToOut(self,translationVector,waveNumber,n)
            
            % Note: translation vector is given as a complex number with
            % real and imaginary parts corresponding to the x and y
            % coordinates respectively
                       
            % process translation vector
            r = abs(translationVector);
            theta=angle(translationVector);

            % in the addition theorem we only sum over some indices so
            % determine those indices here; ii holds the indices that we
            % sum over
            [nn,mm]=meshgrid(-self.nmax:self.nmax,-self.nmax:self.nmax);
            pp=nn-mm;

            % initialise matrix
            matrix = zeros(2*self.nmax+1,2*self.nmax+1);

            % get the necessary Bessel function values
            bess = besselj(0:2*self.nmax,waveNumber*r);
            B=bess(abs(pp)+1);

            % vectorised way to compute the transfer matrix
            matrix = B.*exp(1i*pp*theta).*self.transferMatrix;
            
            % process parameters
            if nargin > 3
                
                if max(n) > self.nmax
                    error('n must not be greater than self.nmax')
                end
                
                % setup indices into matrix
                indexi = (-n(1):n(1)) + self.nmax + 1;
                indexj = (-n(2):n(2)) + self.nmax + 1;
                
                % return only the relevant part of the matrix
                matrix = matrix(indexi,indexj);
                
            end
            
        end
        
        %------------------------------------------------------------------
        % write matrix that effects a rotation
        %------------------------------------------------------------------
        
        function matrix = rotationMatrix(self,angle,n)
            
            % set up diagonal entries on rotation matrix
            tmp = exp(1i*(-self.nmax:self.nmax)*angle);
            
            % create a sparse diagonal matrix with the diagonal entries
            matrix = spdiags(tmp(:),[0],2*self.nmax+1,2*self.nmax+1);
            
            % process parameters
            if nargin > 2
                
                if max(n) > self.nmax
                    error('n must not be greater than self.nmax')
                end
                
                % setup indices into matrix
                indexi = (-n(1):n(1)) + self.nmax + 1;
                indexj = (-n(2):n(2)) + self.nmax + 1;
                
                % return only the relevant part of the matrix
                matrix = matrix(indexi,indexj);
                
            end
            
        end
        
        
    end % end methods
    
end % end classdef