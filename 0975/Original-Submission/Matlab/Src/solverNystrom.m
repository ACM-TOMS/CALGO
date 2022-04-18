% Nystrom sound soft solver
%
%  s = solverNystrom(k,inc,b) creates a solver object with 
%  wavenumber k and incident fields specified in the incident field inc.
%  The scatterer boundary is specified by obstacle object b.
%
%  s = solverNystrom(k,[],b) creates a solver object with wavenumber k and defers
%  setting the incident fields.
%
%  s.setIncidentField(inc) sets the incident field as specified in the cell
%  array inc.
%
%  s.setup(n) sets up the Nystrom discretisation of an integral equation
%  for the sound soft scattering problem using 2n points on the boundary.
%
%  s.solve() solves the Nystrom system for the incident field inc.
%
%  val = s.getFarField(theta) computes the far field of inc{1} at angles
%  specified by theta.
%
%  val = s.getFarField(theta,index) computes the far field of inc{index} at 
%  angles specified by theta. The far field of inc{index(k)} is val(:,k).
%
%  val = s.getCrossSection(self) computes the cross section for the
%  incident field inc{1}.
%
%  val = s.getCrossSection(self,index) computes the cross section for the
%  incident fields inc{index}.
%
%  s.plotFarField(theta) plots the far field for the incident field inc{1}
%  at points with polar angle theta.
%
%  s.plotFarField(theta,index) plots the far field for the incident field 
%  inc{index} at points with polar angle theta.
%
%  s.plotFarField(theta,index,opts) plots the far field for the incident field 
%  inc{index} at points with polar angle theta, with plot linetype
%  specified by opts.
%
%  obj = s.plotFarField(...) returns the handle of the plot.
%
%  s.plotCrossSection(...) plots the cross section in decibels. See
%  plotFarField for the syntax.
%
%  r = s.getRadius() returns the approximate radius r of the obstacle b.
%
%  s.visualize() plots the obstacle b.
%
% See also: solver, solverNystromRobin, mfsSolver.
%
% Stuart C. Hawkins - 20 November 2014

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


classdef solverNystrom < solver
    
    properties
        coupling_parameter
        nystrom_parameter
        matrix
        cof
        scatterer
    end
    
    methods
        
        %-----------------------------------------
        % constructor
        %-----------------------------------------
        
        function self = solverNystrom(kwave,incidentField,scatterer)
            
            %  call parent constructor
            self = self@solver(kwave,incidentField);
            
            % set coupling parameter
            self.coupling_parameter = self.kwave;
            
            % check that the scatterer is of the correct type
            if ~isa(scatterer,'obstacle')
                error('scatterer must be of class obstacle')
            end
            
            % set the scatterer
            self.scatterer = scatterer;
            
            % set some things to empty... this will highlight if
            % setup etc not run
            self.matrix = [];
            self.cof = [];
            
        end
        
        %===============================================================
        % methods defined in the solver class that must be overloaded 
        %===============================================================
        
        %-----------------------------------------
        % setup
        %-----------------------------------------
        
        function setup(self,nystrom_parameter)
            
            % store Nystrom parameter
            self.nystrom_parameter = nystrom_parameter;
            
            % compute the Nystrom matrix
            self.matrix = self.nystrom_matrix(nystrom_parameter);
            
        end
        
        %-----------------------------------------
        % solve
        %-----------------------------------------
        
        function solve(self)
            
            if isempty(self.matrix)
                
                error('Must call setup() first')
                
            end
            
            % setup the right hand side
            b = self.nystrom_rhs(self.nystrom_parameter);
            
            % solve the system
            self.cof = self.matrix \ b;
            
        end
        
        %-----------------------------------------
        % get far field
        %-----------------------------------------
        
        function val = getFarField(self,points,index)
            
            % set default for index
            if nargin<3
                index = 1;
            end
            
            % check that solve has been run
            if isempty(self.cof)
                
                error('Must run solve() first')
                
            end
            
            % compute the far field
            val = self.nystrom_far_field(points,index);
            
            
        end
        
        %===============================================================
        % additional methods that might be useful to the user
        %===============================================================
        
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

            % call the visualize method of the obstacle
            self.scatterer.visualise();
            
        end
        
        %-----------------------------------------
        % get radius
        %-----------------------------------------
        
        function val = getRadius(self)
            
            % evaluate the geometry parametrisation at lots of points and
            % take maximum
            t = 2*pi*(0:999)/1000;
            [x,y,qx,qy] = self.scatterer.geom(t);
            val = sqrt(max(qx.^2+qy.^2));
            
        end
        
    end % end methods
    
    methods(Access=private)
        
        %===============================================================
        % methods that are not for use by the user
        %===============================================================

        %-----------------------------------------
        % setup Nystrom matrix
        %-----------------------------------------
        
        function A = nystrom_matrix(self,n)
            
            % set constant
            C = 0.57721566490153286060651209008240243104215933593992;
            
            %-----------------------------------------------------------
            % set up quadrature for nonsingular integrals
            %-----------------------------------------------------------
            
            pp=pi*(0:2*n-1)/n;
            pp=pp(:);
            pw=pi*ones(2*n,1)/n;
            
            %-----------------------------------------------------------
            % set up quadrature for singular integrals
            %-----------------------------------------------------------
            
            % initialize quadrature matrix
            R=zeros(2*n,2*n);
            
            % setup the quadrature points... see Colton and Kress,
            % Inverse EM and Acoustic Scattering Theory, 3rd Ed, p78
            R(1,:)=-(pi/n^2)*cos(n*(pp(1)-pp));
            
            for m=1:n-1
                R(1,:)=R(1,:) - (2*pi/(n*m))*cos(m*(pp(1)-pp.'));
            end
            
            for i=2:2*n
                
                R(i,:) = circshift(R(1,:),[0,i-1]);
                
            end
            
            %-----------------------------------------------------------
            % set up matrix
            %-----------------------------------------------------------
            
            % initialize matrix
            A = zeros(2*n,2*n);
            
            % apply mapping to quadrature points
            [x,y,qx,qy,dqx,dqy,ddqx,ddqy,nqx,nqy,jac]=self.scatterer.geom(pp);
            
            % Compute the matrix... see Colton and Kress,
            % Inverse EM and Acoustic Scattering Theory, 3rd Ed, p78
            
            for i=1:2*n
                
                d = qx(i)-qx;
                e = qy(i)-qy;
                
                r=sqrt( d.^2 + e.^2 );
                
                % set dummy value for diagonal
                r(i)=1000;
                
                dp = d.*nqx + e.*nqy;
                
                double_kernel = 0.5*1i*self.kwave*besselh(1,self.kwave*r).*dp.*jac./r;
                
                single_kernel = 0.5*1i*besselh(0,self.kwave*r).*jac;
                
                double_sing = -(1/(2*pi))*self.kwave*besselj(1,self.kwave*r).*dp.*jac./r;
                
                single_sing = -(1/(2*pi))*besselj(0,self.kwave*r).*jac;
                
                single_smooth = single_kernel - single_sing.*log(4* (sin(0.5*(pp(i)-pp))).^2 );
                
                double_smooth = double_kernel - double_sing.*log(4* (sin(0.5*(pp(i)-pp))).^2 );
                
                A(i,:) = double_smooth.*pw  + double_sing.*R(:,i) + 1i*self.coupling_parameter*single_smooth.*pw + 1i*self.coupling_parameter*single_sing.*R(:,i);
                
            end
            
            
            % Note the 1/2*pi is cancelled by the quadrature weight and
            % the jacobian is already included
            double_smooth_diagonal=-(dqx.*ddqy-dqy.*ddqx)./(2*n*(dqx.^2+dqy.^2));
            
            single_smooth_diagonal= (0.5i - C/pi - (1/(2*pi))*log( 0.25*self.kwave^2*jac.^2 ) ).*jac.*pw;
            
            single_sing_diagonal = -(1/(2*pi))*besselj(0,0).*jac;
            
            double_sing_diagonal = zeros(size(double_smooth_diagonal));
            
            for i=1:2*n
                A(i,i) = 1 + double_smooth_diagonal(i) + double_sing_diagonal(i)*R(i,i) + 1i*self.coupling_parameter*single_smooth_diagonal(i) + 1i*self.coupling_parameter*single_sing_diagonal(i)*R(i,i);
            end
            
        end
        
        %-----------------------------------------------------------
        % set up right hand side
        %-----------------------------------------------------------
        
        function b = nystrom_rhs(self,n)
            
            % set up quadrature points
            pp=pi*(0:2*n-1)/n;
            pp=pp(:);
            
            % apply mapping to quadrature points
            [x,y,qx,qy]=self.scatterer.geom(pp);
            
            % get qmap in complex format for the incident field evaluation
            qz = qx + 1i*qy;
            
            % setup the right hand side vectors
            for k=1:self.numIncidentField
                
                b(:,k) = -2*self.incidentField{k}.evaluate(qz);
                
            end
            
        end
        
        %-----------------------------------------------------------
        % compute far field
        %-----------------------------------------------------------
        
        function val = nystrom_far_field(self,points,index)
            
            % make sure points is a column vector
            points = points(:);
            
            % initialize temporary array... do this because it is
            % convenient to compute the far field at each point for all
            % right hand sides in one go. This gives the transpose of val
            tmp = zeros(length(index),length(points));
            
            % get coordinates (p,q) of observation point on the unit circle
            p=cos(points);
            q=sin(points);

            % set up quadrature for nonsingular integrals           
            pp=pi*(0:2*self.nystrom_parameter-1)/self.nystrom_parameter;
            pp=pp(:);
            pw=pi*ones(2*self.nystrom_parameter,1)/self.nystrom_parameter;

            % apply parametrisation of boundary
            [x,y,qx,qy,dqx,dqy,ddqx,ddqy,nqx,nqy,jac] = self.scatterer.geom(pp);

            % loop through point
            for j = 1:length(points)
                
                % commpute kernel using details on p75 of Colton and Kress,
                % Inverse EM and Acoustic Scattering Theory, 3rd Ed
                
                d = p(j)-qx;
                e = q(j)-qy;
                
                r=sqrt( d.^2 + e.^2 );
                
                ndp = p(j).*nqx + q(j).*nqy;
                dp = p(j).*qx + q(j).*qy;
                
                kernel = ((1+1i)/sqrt(2))/sqrt(8*pi*self.kwave)*jac.*(-1i*self.kwave*ndp+1i*self.coupling_parameter).*exp(-1i*self.kwave*dp);
                
                % loop through the right hand sides
                for k = 1:length(index)
                
                    % compute the far field for each right hand side
                    tmp(k,j) = sum(pw.*kernel.*self.cof(:,index(k)));
                
                end
                
            end
            
            % take transpose... do this because it is
            % convenient to compute the far field at each point for all
            % right hand sides in one go. This gives the transpose of val
            val = tmp.';

        end
        
    end % end methods
    
end