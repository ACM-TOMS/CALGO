% Circular scatterer
%
%   S = scatterer(x,r,m) creates an object representing a
%   circular scatterer with centre x and radius r. The material properties 
%   of the object are specified by m, which can take the values
%
%                'SOFT' - a Dirichlet BC u = 0 is applied on the boundary
%                'HARD' - a Neumann BC du/dn = 0 is applied on the boundary
%                'ROBIN' - a Robin BC du/dn + 1i * mu * u = 0
%                          is applied on the boundary
%                'DIELECTRIC' - a transmission BC is applied on the 
%                               boundary
%
%   In the case of m = 'ROBIN' the parameter mu must be set eg using the
%   setRobinParameter method. 
%
%   In the case of m = 'DIELECTRIC' the refractive index must be set eg 
%   using the setRefractiveIndex method.
%
%   S = scatterer(x,r,'ROBIN',mu) creates an object representing a
%   circular scatterer with centre x and radius r. A Robin boundary 
%   condition du/dn + 1i * mu * u = 0 is applied on the boundary.
%
%   S = scatterer(x,r,'DIELECTRIC',ri) creates an object representing a
%   circular scatterer with centre x and radius r. A transmission boundary 
%   condition is applied on the boundary and the refractive index in the
%   interior is set to ri.
%
%   S = scatterer(x,r,'DIELECTRIC',ri,rho) creates an object representing a
%   circular scatterer with centre x and radius r. A transmission boundary 
%   condition is applied on the boundary and the refractive index in the
%   interior is set to ri, the density in the interior is set to rho.
%
%   S.setRefractiveIndex(ri) sets the refractive index to ri.
%
%   S.setDensity(rho) sets the density to rho.
%
%   S.setRobinParameter(mu) sets the Robin boundary condition parameter 
%   mu (if the scatterer has a Robin BC).
%
%   ri = S.getRefractiveIndex() gets the refractive index ri. If the
%   scatterer is not dielectric then ri = Inf.
%
%   rho = S.getDensity() gets the density rho. If the
%   scatterer is not dielectric then rho = Inf.
%
%   mu = S.getRobinParameter() gets the Robin boundary condition parameter 
%   mu. If the scatterer does not have a Robin BC then mu = Inf.
%
%   x = S.getCentre() gets the centre of the scatterer.
%
%   r = S.getRadius() gets the radius r of the scatterer.
%
%   s = S.getSize(kwave) gets the acoustic size, which is the diameter of the
%   scatterer divided by the wavelength.
%
%   S.show() plots a visualisation of the scatterer.
%
%   S.show(opts) plots a visualisation of the scatterer with linestyle
%   specified by opts.
%
% Also:
%
%   S.addCoating(r,ri) adds a dielectric layer with refractive index ri to 
%   the scatterer. The outer radius of the layer is r. 
%
%   S.addCoating(r,ri,rho) adds a dielectric layer with refractive index ri 
%   and density rho to the scatterer. The outer radius of the layer is r. 
%
%   S.addCoating may be used several times to create a layered scatterer
%   with a core and with L coatings. This partitions the plane into L+2
%   regions (Region 1 is the core, Region L+2 is the exterior of the
%   scatterer and the other Regions 2,...,L+1 are the layers, numbered from 
%   inside to outside)
%
%   N = S.getNumRegion() returns the number of regions N. Region 1 is the 
%   core and region N is the exterior. Regions 2,...,N-1 correspond to
%   layers.
%
%   S.setRefractiveIndex(j,ri) sets the refractive index of region j to ri.
%
%   S.setDensity(j,rho) sets the density of region j to rho.
%
%   ri = S.getRefractiveIndex() returns a vector ri containing the 
%   refractive index of all regions in the scatterer. (If there are only two
%   regions then only the refractive index of the core is returned.)
%
%   rho = S.getDensity() returns a vector rho containing the density
%   of all regions in the scatterer. (If there are only two
%   regions then only the density of the core is returned.)
%
%   r = S.getRadius(j) gets the outer radius r of Region j.
%
% Note: in the functions above vectors in the plane are represented by
% complex numbers.
%
% See also: MieSolver.
%
% Stuart C. Hawkins - 8 June 2017

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



classdef scatterer < handle
    
    properties (Access=private)
        num_region
        centre
        radius
        refractive_index
        density
        bc
        robin_mu
    end
    
    %=================================================================
    % Public methods
    %=================================================================
    
    methods(Access=public)
        
        %-----------------------------------------------
        % constructor
        %-----------------------------------------------
        
        function self = scatterer(varargin)
            
            % - - - - - - - - - - - - - - - - 
            % copy if nargin==1
            % - - - - - - - - - - - - - - - - 

            if nargin==1
               
                % then varargin{1} is a scatterer and we want to copy it
                self.num_region = varargin{1}.num_region;
                self.centre = varargin{1}.centre;
                self.radius = varargin{1}.radius;
                self.refractive_index = varargin{1}.refractive_index;
                self.density = varargin{1}.density;
                self.bc = varargin{1}.bc;
                self.robin_mu = varargin{1}.robin_mu;
     
                return
                
            end
                
            % pull out the centre, radius and BC from varargin
            centre = varargin{1};
            radius = varargin{2};
            bc = varargin{3};
            
            % - - - - - - - - - - - - - - - - 
            % basic setup
            % - - - - - - - - - - - - - - - - 
 
            % initially a scatterer has only one interior region and the
            % exterior region
            self.num_region=2;
            
            % store the centre and radius
            self.centre=centre;
            self.radius=radius;
            
            % - - - - - - - - - - - - - - - - 
            % process boundary condition
            % - - - - - - - - - - - - - - - - 
            
            if strcmp(upper(bc),'SOFT')
                
                % set the boundary condition
                self.bc=boundary_condition.soft;
                
                % set the corresponding parameters
                self.refractive_index=[Inf 1];
                self.robin_mu=Inf;
                self.density=[Inf 1];
                
            elseif strcmp(upper(bc),'HARD')
                
                % set the boundary condition
                self.bc=boundary_condition.hard;

                % set the corresponding parameters
                self.refractive_index=[Inf 1];
                self.robin_mu=Inf;
                self.density=[Inf 1];

            elseif strcmp(upper(bc),'ROBIN')
                
                % set the boundary condition
                self.bc=boundary_condition.robin;

                % set the corresponding parameters
                self.refractive_index=[Inf 1];
                self.density=[Inf 1];

                if nargin>3
                    % set the Robin BC parameter
                    self.robin_mu = varargin{4};
                else
                    % set Robin BC parameter to NaN... this will need to be
                    % changed later
                    self.robin_mu=NaN;
                end
                
            elseif strcmp(upper(bc),'DIELECTRIC')
                
                % set the boundary condition
                self.bc=boundary_condition.transmission;

                % set the appropriate Robin BC parameter
                self.robin_mu=Inf;

                if nargin>3
                    % set the refractive index
                    self.refractive_index=[varargin{4} 1];
                else
                    % set the refractive index to NaN... this will need to
                    % be changed later
                    self.refractive_index=[NaN 1];
                end
               
                if nargin>4
                    % set the density
                    self.density=[varargin{5} 1];
                else
                    % set the density to NaN... this will need to
                    % be changed later if acoustic transmission BC is used
                    self.density=[NaN 1];
                end
                
            else
                
                error('Material/BC %s not recognised',bc)
                
            end

            % - - - - - - - - - - - - - - - - 
            % display warnings as appropriate
            % - - - - - - - - - - - - - - - - 

            % Note: don't display a warning if the density is NaN... this
            % is only relevant if acoustic transmission conditions are
            % applied and it is checked in mieProblem.addScatterer
            
            % if the refractive index is NaN then display a warning
            if nnz(isnan(self.refractive_index))>0
                warning('Refractive index is not set for some regions. Use setRefractiveIndex method to set.')
            end
                        
            % if the Robin parameter is NaN then display a warning
            if nnz(isnan(self.robin_mu))>0
                warning('Robin parameter is not set for some boundaries. Use setRobinParameter method to set.')
            end
            
        end
        
        %-----------------------------------------------
        % get centre
        %-----------------------------------------------
        
        function val = getCentre(self)
            
            val = self.centre;
            
            
        end

        % and American English version... just a wrapper for the method 
        % getCentre
        function varargout = getCenter(varargin)
            
            varargout{:} = getCentre(varargin{:});
            
        end
        
        %-----------------------------------------------
        % set centre
        %-----------------------------------------------
        
        function setCentre(self,centre)
            
            self.centre = centre;
            
            
        end

        % and American English version... just a wrapper for the method 
        % getCentre
        function setCenter(varargin)
            
            setCentre(varargin{:});
            
        end
        
        %-----------------------------------------------
        % get num region
        %-----------------------------------------------
        
        function val = getNumRegion(self)
            
            val = self.num_region;
            
            
        end
        
        %-----------------------------------------------
        % get refractive index
        %-----------------------------------------------
        
        function val = getRefractiveIndex(self,region)
            
            if nargin<2
                if self.num_region == 2
                    val = self.refractive_index(1);
                else
                    val = self.refractive_index(:);
                end
            else
                val = self.refractive_index(region);
            end
            
        end
        
        %-----------------------------------------------
        % get density
        %-----------------------------------------------
        
        function val = getDensity(self,region)
            
            if nargin<2
                if self.num_region == 2
                    val = self.density(1);
                else
                    val = self.density(:);
                end
            else
                val = self.density(region);
            end
            
        end
        
        %-----------------------------------------------
        % get boundary condition
        %-----------------------------------------------
        
        function val = getbc(self,region)
            
            if nargin==1
                val = self.bc;
            else
                val = self.bc(region);
            end
            
        end
       
        %-----------------------------------------------
        % determine if the scatterer has any transmission
        % boundaries
        %-----------------------------------------------
        
        function val = hasTransmission(self)
            
            val = nnz(self.bc == boundary_condition.transmission) > 0;

        end
        
        %-----------------------------------------------
        % get mu
        %-----------------------------------------------
        
        function val = getRobinParameter(self)
                        
            val = self.robin_mu;
            
        end
                
        %-----------------------------------------------
        % get radius
        %-----------------------------------------------
        
        function rad = getRadius(self,varargin)
            
            % if region is specified then return the radius of that
            % region... otherwise return the radius of the scatterer (which
            % is the radius of the outermost region).
            
            if nargin>1
                region = varargin{1};
            else
                region = length(self.radius);
            end
            
            rad = self.radius(region);
            
            
        end
        
        %-----------------------------------------------
        % get acoustic size
        %-----------------------------------------------
        
        function val = getSize(self,kwave)
            
            % compute wavelength from kwave
            lambda = 2 * pi / kwave;
            
            % divide the diameter by the wavelength
            val = 2 * self.getRadius / lambda;
            
            
        end

        %-----------------------------------------------
        % get mask for a particular region
        %-----------------------------------------------
        
        function val = mask(self,points,region)

            % set padding parameter
            padd = 0;
            
            if region==1
                
                % first region is the disk with radius self.radius(1)
                val = abs(points-self.centre) < (1-padd) * self.radius(1);
                
            elseif region==self.num_region
                
                % last region is outside the disk with radius
                % self.radius(region-1)
                val = abs(points-self.centre) > (1+padd) * self.radius(region-1);

            else
                
                % the region is annular with inner radius
                % self.radius(region-1) and outer radius
                % self.radius(region)
                val = (abs(points-self.centre) < (1-padd) * self.radius(region)) ...
                    & (abs(points-self.centre) > (1+padd) * self.radius(region-1));

            end
               
        end

        %-----------------------------------------------
        % set refractive index
        %-----------------------------------------------

        function self = setRefractiveIndex(self,varargin)
               
            if nargin==3
                % a region has been specified
                region = varargin{1};
                m = varargin{2};
            elseif nargin==2
                % no region has been specified... if there are only two
                % regions then we can assume that the user wants to specify
                % the refractive index for region 1.
                if self.num_region~=2
                    error('The scatterer has %d regions and a region must be specified.',self.num_region);
                end
                region = 1;
                m = varargin{1};
            end                        
            
            % check the specified region is within range
            if region > self.num_region
                
                error('Scatterer does not have region %d.',region)
                
            end
            
            % check that the region is not the outside region
            if region == self.num_region
                
                error('Refractive index for the exterior region may not be set.')
                
            end
            
            % check that the region is dielectric
            if self.bc(region)~=boundary_condition.transmission
            
                error('Region %d is not dielectric.',region)

            end
                
            % set refractive index
            self.refractive_index(region) = m;
            
        end
        
        %-----------------------------------------------
        % set density
        %-----------------------------------------------

        function self = setDensity(self,varargin)
               
            if nargin==3
                % a region has been specified
                region = varargin{1};
                m = varargin{2};
            elseif nargin==2
                % no region has been specified... if there are only two
                % regions then we can assume that the user wants to specify
                % the density for region 1.
                if self.num_region~=2
                    error('The scatterer has %d regions and a region must be specified.',self.num_region);
                end
                region = 1;
                m = varargin{1};
            end                        
            
            % check the specified region is within range
            if region > self.num_region
                
                error('Scatterer does not have region %d.',region)
                
            end
            
            % check that the region is not the outside region
            if region == self.num_region
                
                error('Density for the exterior region may not be set.')
                
            end
            
            % check that the region is dielectric
            if self.bc(region)~=boundary_condition.transmission
            
                error('Region %d is not dielectric.',region)

            end
                
            % set density
            self.density(region) = m;
            
        end
        
        %-----------------------------------------------
        % set robin BC parameter
        %-----------------------------------------------

        function self = setRobinParameter(self,mu)
            
            % check that the scatterer is specified to have a Robin BC
            if self.bc(1) ~= boundary_condition.robin
                
                error('scatterer does not have Robin BC.')
                
            end
            
            % set the Robin BC parameter
            self.robin_mu = mu;
            
        end
        
        %-----------------------------------------------
        % add a region
        %-----------------------------------------------
        
        % The regions are indexed starting from the centre and the last
        % region is the exterior region.
        %
        % When a region is added, it must be inserted between the previous
        % outer shell and the exterior region.
        
        function self = addCoating(self,outer_radius,refractive_index,density)
        
            % set default for density
            if nargin<4
                % NaN means density is not set... this parameter has
                % meaning only for acoustic transmission BCs
                density = NaN;
            end
            
            % check that the outer_radius is outside the rest of the 
            % scatterer
            if outer_radius <= max(self.radius)
                
                error('outer_radius must be greater than max(self.radius).')
                
            end
            
            % check that the refractive index is finite
            if ~isfinite(refractive_index)
                
                error('refractive_index of coating must be finite.')
                
            end
            
            % update scatterer details
            self.num_region = self.num_region+1;

            % move exterior region details outwards
            self.refractive_index(self.num_region) = self.refractive_index(self.num_region-1);
            self.density(self.num_region) = self.density(self.num_region-1);
            
            % add new region details
            self.refractive_index(self.num_region-1) = refractive_index;
            self.density(self.num_region-1) = density;
            self.radius(self.num_region-1) = outer_radius;
            self.bc(self.num_region-1) = boundary_condition.transmission;
            
        end
        
        %-----------------------------------------------
        % overload disp
        %-----------------------------------------------
            
        function disp(self)
           
            fprintf('Circular scatterer object center (%4.2f,%4.2f)\n\n',real(self.centre),imag(self.centre));
            fprintf('Regions\n');
            fprintf('-------\n');
            if self.bc(1)==boundary_condition.soft
                str='[BC SOFT]';
            elseif self.bc(1)==boundary_condition.hard
                str='[BC HARD]';
            elseif self.bc(1)==boundary_condition.robin
                str=sprintf('[BC ROBIN (%4.3e)]',self.robin_mu);
            else
                str='';
            end
            fprintf('(%4.2f,%4.2f) ri=%4.3f density=%4.3f %s\n',0,self.radius(1),self.refractive_index(1),self.density(1),str);
            for k=2:self.num_region-1
                fprintf('(%4.2f,%4.2f) ri=%4.3f density=%4.3f \n',self.radius(k-1),self.radius(k),self.refractive_index(k),self.density(k));
            end
            fprintf('(%4.2f,%4.2f) ri=%4.3f density=%4.3f \n',self.radius(self.num_region-1),Inf,self.refractive_index(self.num_region),self.density(self.num_region));
                       
        end
        
        %-----------------------------------------------
        % show scatterer
        %-----------------------------------------------

        % Plot the scatterer, with some basic information about its
        % material properties.
        
        function obj_out = show(self,opts,label)
                                  
            % set default for opts
            if nargin<2 || isempty(opts)
                opts = 'k-';
            end
            
            % set default for label
            if nargin<3
                label = 1;
            end
            
            % initialise graphic object
            obj=[];
            
            % create array of angles
            theta=(0:100)*2*pi/100;
            
            % create points on unit circle
            z=exp(1i*theta);
           
            % store hold state and put hold on
            store_hold = ishold;
            hold on
            
            % plot the circlular boundaries between regions
            for k=1:self.num_region-1
                
                % compute points on the circle
                w = self.centre + self.radius(k)*z;
                
                % plot the points
                tobj=plot(real(w),imag(w),opts);
                
                % store the graphics handle to the new line
                obj=[obj;tobj];

            end
            
            % add labels if required
            if label~=0
            
                % add labels for refractive indices
                for k=1:self.num_region-1
                    
                    % compute coordinates of the centre of the region
                    if k==1
                        r=self.centre+0.5*self.radius(k);
                    else
                        r=self.centre+0.5*(self.radius(k)+self.radius(k-1));
                    end
                    
                    % get string to print
                    str=sprintf('ri=%3.2f',self.refractive_index(k));
                    
                    % print the string
                    tobj=text(real(r),imag(r),str);
                    set(tobj,'backgroundcolor','white');
                    
                    % store the graphics handle to the text
                    obj=[obj;tobj];
                    
                end
                
            end
            
            % restore hold state
            if ~store_hold
                hold off
            end
            
            % setup return object if needed
            if nargout>0
                obj_out = obj;
            end
            
        end
    
    end
    
end