% Solver
%
%  s = solver(k,inc) creates a solver object with wavenumber k and incident
%  fields specified in the incident field inc.
%
%  s = solver(k,[]) creates a solver object with wavenumber k and defers
%  setting the incident fields.
%
%  s.setIncidentField(inc) sets the incident field as specified in the cell
%  array inc.
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
% Warning: this should be considered an ABSTRACT base class... it is not
% intended that objects of this class be instantiated. Several methods of
% this class must be overridden.
%
% See also: solverNystrom, solverNystromRobin, mfsSolver.
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


classdef solver < handle
    
    properties
        kwave
        incidentField
        numIncidentField
    end
    
    methods
        
        %===============================================================
        % methods that will (probably) not be overridden in child class
        %===============================================================

        %-----------------------------------------
        % constructor
        %-----------------------------------------
        
        function self = solver(kwave,incidentField)
            
            self.kwave = kwave;
            
            if ~iscell(incidentField)
                
                self.incidentField{1} = incidentField;
               
            else
                
                self.incidentField = incidentField;
                
            end
            
            self.numIncidentField = length(self.incidentField);
            
        end
        
        %-----------------------------------------
        % set incident field
        %-----------------------------------------
        
        function setIncidentField(self,incidentField)
            
            % make sure that the incident field is stored as a cell array
            % even if there is only one incident field
            if ~iscell(incidentField)
                
                self.incidentField{1} = incidentField;
               
            else
                
                self.incidentField = incidentField;
                
            end
            
            % set number of incident fields
            self.numIncidentField = length(self.incidentField);
            
        end
        
        %-----------------------------------------
        % get cross section
        %-----------------------------------------

        function val = getCrossSection(self,points,index)
            
            % set default for index
            if nargin < 3
                index = 1;
            end
            
            val = 10*log10(2*pi*abs(self.getFarField(points,index)).^2);
            
        end
        
        %-----------------------------------------
        % plot far field
        %-----------------------------------------
        
        function varargout = plotFarField(self,points,index,opts)
                        
            % set default for opts
            if nargin<4
                opts = 'b-';
            end
            
            % set default index
            if nargin < 3 || isempty(index)
                index = 1;
            end
            
            % get the farfield
            val = self.getFarField(points,index);
            
            % do the plot
            obj = plot(points,real(val),opts,points,imag(val),opts);
            
            % return value if necessary
            if nargout > 0
                varargout{1} = obj;
            end
            
        end
        
        %-----------------------------------------
        % plot cross section
        %-----------------------------------------
        
        function varargout = plotCrossSection(self,points,index,opts)
                        
            % set default for opts
            if nargin<4
                opts = 'b-';
            end
            
            % set default index
            if nargin < 3 || isempty(index)
                index = 1;
            end
            
            % get the farfield
            val = self.getFarField(points,index);
            
            % do the plot
            obj = plot(points,10*log10(2*pi*abs(val).^2),opts);
            
            % return value if necessary
            if nargout > 0
                varargout{1} = obj;
            end
            
        end
        
    end % end methods
    
    methods(Abstract=true)
        
        %===============================================================
        % these methods must be overridden in the child class
        %===============================================================

        %-----------------------------------------
        % setup
        %-----------------------------------------
        
        setup(self)
        
        %-----------------------------------------
        % solve
        %-----------------------------------------
        
        solve(self)
        
        %-----------------------------------------
        % get far field
        %-----------------------------------------
        
        val = getFarField(self,points,index)
        
    end % end methods
    
end