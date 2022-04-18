% Cylindrical Mie Series solver
%
%   P = MieSolver(inc) where inc is of class 'incident' creates an
%   instance of the Cylindrical Mie Series solver.
%
%   P = MieSolver(inc,S1,...,SN) where inc is of class 'incident' creates
%   an instance of the Cylindrical Mie Series solver with scatterers
%   S1,...,SN.
%
%   P = MieSolver(inc,{S1,...,SN}) where inc is of class 'incident'
%   creates an instance of the Cylindrical Mie Series solver with 
%   scatterers S1,...,SN.
%
%   P.addScatterer(S) where S is of class 'scatterer' adds the scatterer S
%   to the scattering configuration.
%
%   P.transmissionTM() configures the solver so that any transmission
%   boundary conditions occurring are of TM type (field is continuous and
%   the normal derivative divided by refractive index squared are
%   continuous).
%
%   P.transmissionTE() configures the solver so that any transmission
%   boundary conditions occurring are of TE type (field and its normal
%   derivative are continuous).
%
%   P.transmissionAcoustic() configures the solver so that any transmission
%   boundary conditions occurring are of acoustic type (field is continuous
%   and its normal derivative divided by density are continuous).
%
%   P.transmissionCustom(f) configures the solver so that the transmission
%   boundary conditions applies f(m) * normal derivative of the field is
%   continuous across the interface. Here m is the refractive index.
%
%   S = P.getScatterer(j) returns the jth scatterer in the scattering
%   configuration.
%
%   P.solve() computes the Cylindrical Mie Series coefficients for the
%   scattered wave. The default truncation parameter for the series uses
%   Wiscombe's formula.
%
%   P.incrementNmax(d) increments the truncation parameter for the series.
%
%   n = P.getNmax() returns the vector of truncation parameters for the 
%   series. The jth entry is the truncation parameter for the jth
%   scatterer.
%
%   n = P.getNmax(j) returns the truncation parameter for jth scatterer.
%
%   u = P.getInducedField(x) returns the induced field u at the points
%   specified by x. For points x in the exterior region, u is the scattered
%   field. For points x in interior dielectric regions, u is the total
%   field. For points x in interior inpenetrable regions (ie the core with
%   SOFT, HARD or ROBIN boundary conditions), u is zero.
%
%   u = P.getTotalField(x) returns the total field u at the points
%   specified by x. For points x in interior inpenetrable regions (ie the 
%   core with SOFT, HARD or ROBIN boundary conditions), u is zero. In other
%   regions u is the total field.
%
%   u = P.getFarfield(x) returns the far field u at the points
%   specified by x. It is assumed that the points in x lie on the unit
%   circle ie that abs(x) = ones(size(x)).
%
%   u = P.getRcs(x) returns the cross section u at the points
%   specified by x. It is assumed that the points in x lie on the unit
%   circle ie that abs(x) = ones(size(x)).
%
%   P.visualiseTotalField() visualises the real part of the total field in
%   a rectangular region containing the scatterers with 10 points per
%   (exterior) wavelength.
%
%   P.visualiseTotalField(lim) visualises the real part of the total field 
%   in a rectangular region with x in [min(real(lim)),max(real(lim))] and
%   y in [min(imag(lim)) max(imag(lim))].
%
%   P.visualiseTotalField(lim,opts) visualises the total field using
%   options specified in opts. When opts includes 'abs' the absolute value
%   of the total field is plotted. When opts includes '-f' the limt of
%   400^2 points in the grid is over-ruled.
%
%   P.visualiseRcs() visualises the cross section of the configuration.
%
%   P.visualiseRcs(opts) visualises the cross section of the configuration
%   with linestype specified in opts.
%
% Note: in the functions above vectors in the plane are represented by
% complex numbers.
%
% See also: scatterer.
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


classdef MieSolver < handle
    
    properties(Access=public)
        num_scatterers
        scatterer
        nmax
        incident_field
        kwave
        cof
        transmissionType
        transmissionTypeHasBeenSet
        transmissionFun
    end
    
    %=================================================================
    % methods with standard access
    %=================================================================
    
    methods

        %-----------------------------------
        % set the kind of transmission BC
        % used
        %-----------------------------------

        function set.transmissionType(self,type)
            
            switch type
                
                case 'TE'
                    
                    self.transmissionType = 'TE';
                    
                case 'TM'
                    
                    self.transmissionType = 'TM';
                    
                case 'acoustic'
                    
                    self.transmissionType = 'acoustic';
                    
                case 'custom'
                    
                    self.transmissionType = 'custom';
                                        
                otherwise
                    
                    error('Transmission type %s not recognised',type)
                    
            end
           
            % record that the user has set the transmission type
            self.transmissionTypeHasBeenSet = 1;
            
        end

    end
    
    %=================================================================
    % Public methods
    %=================================================================
    
    methods(Access=public)
        
        %-----------------------------------
        % Constructor. Sets up a problem object
        % given an incident field. The wavenumber 
        % is picked up from the incident field
        %-----------------------------------

        function self = MieSolver(incident_field,varargin)
            
            % set the incident field
            self.incident_field = incident_field;
            
            % get the wavenumber from the incident field
            self.kwave = incident_field.kwave;
            
            % set the other things to empty/zero to show no scatterers are
            % added and no coefficients have been calculated
            self.scatterer={};
            self.num_scatterers = 0;
            self.cof=[];
            self.nmax=[];
            
            % set default transmission BC type to TE
            self.transmissionType = 'TE';
            
            % set default transmission function to empty
            self.transmissionFun = [];
            
            % record that transmission type has not been manually set...
            % this allows a warning to be displayed if a DIELECTRIC
            % scatterer is added
            self.transmissionTypeHasBeenSet = 0;
            
            % if other parameters are given, they are scatterers and should
            % be added....
            if nargin>1
                
                % put the scatterers into a list
                if iscell(varargin{1})
                    
                    list = varargin{1};
                
                    % if varargin{1} is a cell then we will ignore later
                    % inputs... display a warning
                    if nargin>2
                        warning('Some inputs were ignored.')
                    end
                    
                else
                    
                    list = {varargin{:}};
                    
                end
                    
                for j = 1:length(list)
               
                    % use the addScatterer method, which performs various
                    % checks on the scatterers, including whether they
                    % intersect
                    self.addScatterer(list{j});
                    
                end
                
            end
            
        end
       
        %-----------------------------------
        % increment nmax
        %-----------------------------------

        % Default is to use Wiscombe's formula for the series truncation
        % parameter nmax. Sometimes more (or less) accuracy is needed and
        % this can be obtained by incrementing nmax.
        
        function incrementNmax(self,increment)
            
            % ensure nmax is set before using this
            if isempty(self.nmax)
                
                error('nmax is not set. This is probably because no scatterers have been added.')
                
            end
            
            % increment nmax
            self.nmax = self.nmax + increment;
            
        end
            
        %-----------------------------------
        % get nmax
        %-----------------------------------

        function val = getNmax(self,scatterer)
        
            if nargin==1 || isempty(scatterer) 
                val = self.nmax;
            else
                val = self.nmax(scatterer);
            end
            
        end
                
        %-----------------------------------
        % The problem object initially has no
        % scatterers. These are added manually using
        % this method.
        %-----------------------------------

        function varargout = addScatterer(self,scatterer)

            % check scatterer is of the correct class
            if ~isa(scatterer,'scatterer')
                error('Scatterer must be of class scatterer.')
            end
            
            % check the scatterer has appropriate value for refractive
            % index... any NaNs suggest material is DIELECTRIC but the
            % refractive index is not set
            if nnz(isnan(scatterer.getRefractiveIndex))>0
                error('Scatterer''s refractive index is not set for some regions. Use setRefractiveIndex method to set.')
            end
            
            % check the scatterer has appropriate value for Robin
            % parameter... any NaNs suggest boundary has Robin BC but the
            % parameter is not set
            if nnz(isnan(scatterer.getRobinParameter))>0
                error('Scatterer''s Robin parameter is not set for some boundaries. Use setRobinLambda method to set.')
            end
            
            % check if the scatterer has any transmission BCs... if it does
            % and if the transmission type has not been set then display a
            % warning
            if scatterer.hasTransmission() && ~self.transmissionTypeHasBeenSet 
                warning('Using default TE type transmission BCs for DIELECTRIC regions')                
            end
            
            % check that, if the scatterer has transmission BCs and the
            % transmission type is acoustic, then the density has been set
            if scatterer.hasTransmission() && strcmp(self.transmissionType,'acoustic')...
              && nnz(isnan(scatterer.getDensity))>0
              
              error(['Using acoustic type transmission BCs for DIELECTRIC regions',...
                  ' but the density has not been set for this scatterer.'])

            end
            
            % check distance between this scatterer and the
            % other scatterers
            for j = 1: self.num_scatterers
            
                % compte the distance
                d =  abs(scatterer.getCentre() - self.scatterer{j}.getCentre()) ...
                    - scatterer.getRadius() - self.scatterer{j}.getRadius();
                
                % check for distance < 0 which indicates scatterers overlap
                if d<0
                    error('Scatterer overlaps existing scatterer %d.',j)
                end
                
            end
            
            % store the scatterer
            self.num_scatterers = self.num_scatterers + 1;
            self.scatterer{self.num_scatterers}=scatterer;
            
            % compute nmax for this scatterer based on the radius of the scatterer
            self.nmax(self.num_scatterers)=wismieval(scatterer.getSize(self.kwave));
            
            % set return value to scatterer index if necessary
            if nargout > 0
                varargout{1} = self.num_scatterers;
            end
            
        end
        
        %-----------------------------------
        % set the transmission BC type
        %-----------------------------------

        function transmissionTE(self)
            
            self.transmissionType = 'TE';
            
        end
        
        %-----------------------------------
        % set the transmission BC type
        %-----------------------------------

        function transmissionTM(self)
            
            self.transmissionType = 'TM';
            
        end
        
        %-----------------------------------
        % set the transmission BC type
        %-----------------------------------

        function transmissionAcoustic(self)
            
            self.transmissionType = 'acoustic';
            
        end
        
        %-----------------------------------
        % set the transmission BC type
        %-----------------------------------

        function transmissionCustom(self,fun)
            
            self.transmissionType = 'custom';
            self.transmissionFun = fun;
            
        end
        
        %-----------------------------------
        % This routine must be called before
        % the exterior field etc are called.
        % This is the routine that sets up
        % and solves the matrix equation.
        %-----------------------------------

        function solve(self)
           
            % compute induced field coefficients
            self.cof = self.solveSystem();
                        
        end

        %-----------------------------------
        % Get the k-th scatterer
        %-----------------------------------

        function val = getScatterer(self,k)
            
            val = self.scatterer{k};
            
        end

            
        %-----------------------------------
        % Get the field induced by a particular 
        % scatterer at points.
        %-----------------------------------

        function val = getInducedField(self,points)
                        
            % check that the problem has been solved
            if isempty(self.cof)
                error('Must call problem.solve first')
            end
            
            % initialise return aray
            val=zeros(size(points));
            
            % get index data for the matrix/coefficient vector
            [nm,nb,aa] = self.matrixDimensions();
            
            % - - - - - - - - - - - - - - - - - - - - -
            % get induced field in the interior regions
            % - - - - - - - - - - - - - - - - - - - - -
                        
            % loop through the scatterers
            for scatk = 1:self.num_scatterers
                
                % loop through the regions of the scatterer
                for k=1:self.scatterer{scatk}.getNumRegion()-1
                    
                    % generate the mask for the region
                    mask = self.scatterer{scatk}.mask(points,k);
                    
                    % compute the induced field on the interior masked
                    % regions
                    val = val + self.get_induced_field_region( ...
                        scatk,k,self.cof(aa{scatk}),points,mask);
                                       
                end
                
            end
            
            % - - - - - - - - - - - - - - - - - - - - -
            % get induced field in the exterior region
            % - - - - - - - - - - - - - - - - - - - - -
            
            % initialise mask
            mask = ones(size(points));
            
            % loop through the scatterers to assemble the mask
            for scatk = 1:self.num_scatterers
                
                % the exterior region is always the last region
                k=self.scatterer{scatk}.getNumRegion();
                    
                % the mask is the intersection of the exterior of all the masks
                mask = mask & self.scatterer{scatk}.mask(points,k);
                    
            end
                
            % loop through the scatterers
            for scatk = 1:self.num_scatterers

                % the exterior region is always the last region
                k=self.scatterer{scatk}.getNumRegion();
                
                % compute the induced field from scatterer scatk in the
                % exterior region
                val = val + self.get_induced_field_region(...
                    scatk,k,self.cof(aa{scatk}),points,mask);
                
            end
            
        end
        
        %-----------------------------------
        % Get the incident field
        %-----------------------------------

        function val = getIncidentField(self,points)

            % get the mask for the exterior region so that we can
            % add in the incident field
            % initialise mask
            mask = ones(size(points)); 
            for scatk = 1:self.num_scatterers                
                k = self.scatterer{scatk}.getNumRegion;
                mask = mask & self.scatterer{scatk}.mask(points,k);
            end
            
            val = self.incident_field.evaluate(points,mask);
            
        end
            
        %-----------------------------------
        % Get the total field induced by
        % the scatterer. This is the induced
        % field in the interior regions but 
        % is the induced field plus the incident
        % field in the exterior.
        %-----------------------------------

        function val = getTotalField(self,points)
            
            % check that the problem has been solved
            if isempty(self.cof)
                error('Must call problem.solve first')
            end
            
            % get the induced field
            val = self.getInducedField(points);
                        
            % compute the incident field and add it on
            val = val + self.getIncidentField(points);
            
        end
                
        %-----------------------------------
        % Get the farfield induced by the scatterer.
        %-----------------------------------

        function val = getFarfield(self,points)
            
            % check that the problem has been solved
            if isempty(self.cof)
                error('Must call problem.solve first')
            end
            
            % get index data for the matrix/coefficient vector
            [nm,nb,aa] = self.matrixDimensions();
            
            % initialise return array
            val=zeros(size(points));
            
            % loop through scatterers
            for scatk = 1:self.num_scatterers
            
                % get the key to the different variables
                % (this lets us know which correspond to the exterior field
                % and interior fields etc)
                variables = self.get_variables(self.scatterer{scatk});

                % select the last region - this is the exterior region
                region=self.scatterer{scatk}.getNumRegion;
                if variables(2*region)=='H'
                    
                    % compute the coefficients corresponding to the
                    % exterior region
                    ii=((2*region-1)*nb(scatk)+1):((2*region)*nb(scatk));

                    % compute the far field
                    val = val + sumcof(points,self.scatterer{scatk}.getCentre(),...
                        self.kwave,self.cof(aa{scatk}(ii)),'F');
                end
                
            end
                        
        end
        
        %-----------------------------------
        % Get the RCS induced by the scatterer.
        %-----------------------------------

        function val = getRcs(self,points)
            
            % check that the problem has been solved
            if isempty(self.cof)
                error('Must call problem.solve first')
            end
            
            % compute the RCS frim the farfield
            val = 10*log10(2*pi*abs(self.getFarfield(points)).^2);
            
        end
        
        %-----------------------------------
        % visualise the total field
        %-----------------------------------

        % wrapper for US spelling
        
        function varargout = visualizeTotalField(varargin)
            
            [varargout{1:nargout}] = visualiseTotalField(varargin{:});

        end
        
        %-----------------------------------
        % visualise the total field
        %-----------------------------------

        function visualiseTotalField(self,lim,opts)
           
            % set default for opts
            if nargin<3                
                opts = '';
            end
            
            % set default lim based on position of scatterers
            if nargin<2 || isempty(lim)
                
                % initialise arrays
                centres = zeros(self.num_scatterers,1);
                radii = zeros(self.num_scatterers,1);
                
                % get the centres and radii of the scatterers
                for j = 1:self.num_scatterers
                    centres(j) = self.scatterer{j}.getCentre();
                    radii(j) = self.scatterer{j}.getRadius();
                end
                
                % set lim array to go from the bottom-left-most scatterer
                % to the top-right-most scatterer           
                lim = [
                    min(real(centres))-max(radii) + 1i*(min(imag(centres))-max(radii))
                    max(real(centres))+max(radii) + 1i*(max(imag(centres))+max(radii))
                    ];                    
                
                % enlarge lim by 50% about its centere
                lim = mean(lim) + 1.5*(lim - mean(lim));
                
            end
            
            % compute with of plot region in both directions
            du = max(real(lim))-min(real(lim));
            dv = max(imag(lim))-min(imag(lim));
            
            % get wavelength... will measure width of plot region in
            % wavelelengths
            lambda = self.wavelength();
            
            % work out how many points are needed in each direction based
            % on the wavelength
            nu = max(100,ceil(10 * du / lambda));
            nv = max(100,ceil(10 * dv / lambda));
            
            % check if number of points is large
            if (nu*nv > 400^2) && ~contains(opts,'-f')
                
                error(['Total number of points is %d > 400^2 and calculation '...
                    'may take a long time. Run with opts=''-f'' to force.'],nu*nv);
                
            end
                
            % set up grid... u will be used for x values and v will be used
            % for y values
            u = linspace(min(real(lim)),max(real(lim)),nu);
            v = linspace(min(imag(lim)),max(imag(lim)),nu);
            
            % generate the 2D grid
            [x,y] = meshgrid(u,v);            
            
            % represent points in the grid in complex format
            z = x + 1i * y;
            
            % get total field
            u = self.getTotalField(z);

            % plot the field
            if contains(opts,'abs')
                % plot the absolute value
                surf(x,y,zeros(size(u)),abs(u))
            else
                % plot the real part
                surf(x,y,zeros(size(u)),real(u))
            end
            
            % store the hold state
            hold_state = ishold;
            
            % add the scatterers
            for j=1:self.num_scatterers

                % plot the scatterer
                hold on
                self.scatterer{j}.show()
                                
            end
            
            % restore hold state if necessary
            if ~hold_state
                hold off
            end
            
        end
        
        %-----------------------------------
        % visualise the RCS
        %-----------------------------------

        % wrapper for US spelling
        
        function varargout = visualizeRcs(varargin)
            
            [varargoute] = visualiseRcs(varargin{:});

        end
        
        %-----------------------------------
        % visualise the RCS
        %-----------------------------------

        function varargout = visualiseRcs(self,opts)
        
            % set default for opts
            if nargin<2
                opts = 'k-';
            end
            
            % get points on which to get the RCS
            theta = linspace(0,2*pi,1000);
            
            % compute the RCS
            u = self.getRcs(exp(1i*theta));
            
            % plot
            obj = plot(theta,u,opts);
            
            % return graphics handle if required
            if nargout > 0
                varargout{1} = obj;
            end           
            
        end
        
        %-----------------------------------
        % visualise the Mie matrix
        %-----------------------------------

        function spyMieMatrix(self)
            
            % get the Mie matrix
            [matrix,indexes,localIndexes,rows,columns,bcrows] = self.setupMieMatrix();
            
            % get dimensions of the various blocks in the Mie matrix
            % nm - blocks corresponding to each scatterer
            % nm - blocks corresponding to regions in each scatterer
            [nm,nb,aa] = self.matrixDimensions();
            
            % compute the indexes corresponding to the end of the blocks
            % corresponding to each scatterer
            snm = cumsum(nm);
            
            % compute the indexed corresponding to the end of the blocks
            % corresponding to each region in the first scatterer
            iblock = (nb(1):nb(1):nm(1)).';

            for j=2:length(nm)
                
                % compute the indexed corresponding to the end of the blocks
                % corresponding to each region in the jth scatterer
                iblock = [iblock;(snm(j-1)+nb(j):nb(j):snm(j)).'];
                    
            end

            % Mie matrix is setup so that its block structure is same in
            % the other dimension ie the blocks are all square
            jblock = iblock;
                        
            % use spyb to visualise the nonzeros and the block structure
            spyb(matrix,iblock,jblock)
            
        end
        
        %-----------------------------------
        % visualise the Multiple scattering
        % matrix
        %-----------------------------------

        function spyMultipleMatrix(self)
          
            % get the Mie matrix
            [matrix,indexes,localIndexes,rows,columns,bcrows] = self.setupMieMatrix();
            
            % get the multiple scattering matrix
            T = self.setupMultipleMatrix(localIndexes);
            
            % get dimensions of the various blocks in the Mie matrix
            % nm - blocks corresponding to each scatterer
            % nm - blocks corresponding to regions in each scatterer
            [nm,nb,aa] = self.matrixDimensions();
            
            % compute the indexes corresponding to the end of the blocks
            % corresponding to each scatterer
            snm = cumsum(nm);
            
            % compute the indexed corresponding to the end of the blocks
            % corresponding to each region in the first scatterer
            iblock = (nb(1):nb(1):nm(1)).';

            for j=2:length(nm)
                
                % compute the indexed corresponding to the end of the blocks
                % corresponding to each region in the jth scatterer
                iblock = [iblock;(snm(j-1)+nb(j):nb(j):snm(j)).'];
                    
            end

            % Mie matrix is setup so that its block structure is same in
            % the other dimension ie the blocks are all square
            jblock = iblock;
                       
            % use spyb to visualise the nonzeros and the block structure
            spyb(T,iblock,jblock)
            
        end
        
    end
    
    %=================================================================
    % Private methods
    %=================================================================
    
    methods %(Access=protected)
        
        %-----------------------------------
        % get wavelength based on wavenumber
        %-----------------------------------

        function val = wavelength(self)
            
            val = 2*pi/self.kwave;
            
        end
        
        %-----------------------------------
        % internal routine to setup
        % the Mie system matrix
        %-----------------------------------
        
        % This is the function signature, to ensure it is protected. The
        % function definition is in setupMieMatrix.m
        
        [matrix,indexes,localIndexes,rows,columns,bcrows] = setupMieMatrix(self);
               
        %-----------------------------------
        % internal routine to setup multiple
        % scattering transformation
        %-----------------------------------

        function matrix = setupMultipleMatrix(self,localIndexes)
            
            % create translation addition theorem object
            taobj = translation(max(self.nmax));

            % work out dimensions of matrix and setup the matrix
            % in sparse format
            [nm,nb,aa] = self.matrixDimensions();
            matrix = speye(sum(nm));

            % setup indexes
            for scatk = 1:self.num_scatterers
                for j=1:2*self.scatterer{scatk}.getNumRegion()
                    localIndexes{scatk}{j} = (j-1)*nb(scatk)+1:j*nb(scatk);
                end
            end
            
            % insert the translation-addition theorem part
            for scatj = 1:self.num_scatterers
                for scati = 1:self.num_scatterers
                    
                    if scatj~=scati
                        
                        % compute the translation vector
                        v = self.scatterer{scati}.getCentre() - ...
                            self.scatterer{scatj}.getCentre();
                        
                        % compute the translation matrix
                        matrix(aa{scati}(localIndexes{scati}{end-1}),aa{scatj}(localIndexes{scatj}{end})) = ...
                             taobj.translationMatrixOutToIn(v,self.kwave,[self.nmax(scati),self.nmax(scatj)]);
                        
                    end
                    
                end
            end
            
        end
        
        %-----------------------------------
        % internal routine to setup and solve
        % the Mie system
        %-----------------------------------
        
        function cof = solveSystem(self)
                        
            % compute Mie matrix with values of the various wavefunctions
            % at the interfaces, also various indexed
            [A,indexes,localIndexes,rows,columns,bcrows] = self.setupMieMatrix();
            
            % if there are more than a single scatterer then we need to
            % compute the multiple scattering matrix
            if self.num_scatterers > 1
                
                % set up multiple scattering matrix
                T = self.setupMultipleMatrix(localIndexes);
                
                % apply the transformation
                A = A*T;
                
            end
            
            % setup incident field coefficients
            tmp = cumsum([0,2*self.nmax+1]);
            inc_cof = zeros(tmp(end),1);
            for scatk=1:self.num_scatterers
                inc_cof(tmp(scatk)+1:tmp(scatk+1)) = ...
                    self.incident_field.get_coefficients(self.scatterer{scatk}.getCentre(),self.nmax(scatk));
            end

            % assemble the reduced matrix and RHS
            L=A(rows,columns);
            b=-A(rows,bcrows)*inc_cof;                       
            
            % solve the linear system.
            % Note that the other parts of cof are set to zero by default.
            % This corresponds to zero field arising from unused
            % coefficients eg the H part in the central region has zero
            % coefficients because the H basis is inadmissable in the
            % central region.
            cof=zeros(size(A,1),1);
            cof(columns)=L \ b;
            
        end
        
        %-----------------------------------
        % Produces a key to the meaning of the
        % unknowns in cof. Each region has two entries
        % A J indicates coefficients of Bessel fns.
        % A H indicates coefficients of Hankel fns.
        % A 0 indicated these variables are set to 0.
        % An I indicates these variables are specified by the
        % incident field.
        %-----------------------------------

        function variables = get_variables(self,scatterer)
            
            % intialise all regions to have J and H
            variables=repmat(['J';'H'],scatterer.getNumRegion(),1);
            
            % first region H part is set to 0 because H unbounded
            % at 0 in the first region
            variables(2)='0';
            
            % J part of the last region is set to I because these are
            % set by the incident field
            variables(end-1)='I';
            
            % if the first region has infinite refractive index
           % then it has J part set to zero as well.
            if ~isfinite(scatterer.getRefractiveIndex(1))
                variables(1)='0';
            end
            
        end
        
        %-----------------------------------
        % Get induced field in a region
        % at points specified by points and mask.
        %-----------------------------------

        function val = get_induced_field_region(self,scattererIndex,region,cof,points,mask)
            
            % initalise return value
            val=zeros(size(points));
            
            % if a mask is specified then apply the mask
            if nargin>5
                points=points(mask);
            end
            
            % check for empty region
            if isempty(points)
                return
            end
            
            % get the key to the variables in cof
            variables=self.get_variables(self.scatterer{scattererIndex});
            
            % initialise a dummy return value. We use a dummy so that
            % the mask can be applied and the computed values inserted into
            % val
            v=zeros(size(points));
            
            % get the size of a block
            nb=2*self.nmax(scattererIndex)+1;
            
            % compute the sum of the J part in the specified region            
            if variables(2*region-1)=='J'
                ii=((2*region-2)*nb+1):((2*region-1)*nb);
                kwave=self.kwave*self.scatterer{scattererIndex}.getRefractiveIndex(region);
                v = v + sumcof(points,self.scatterer{scattererIndex}.getCentre,kwave,cof(ii),'J');
            end
            
            % compute the sum of the H part in the specified region            
            if variables(2*region)=='H'
                ii=((2*region-1)*nb+1):((2*region)*nb);
                kwave=self.kwave*self.scatterer{scattererIndex}.getRefractiveIndex(region);
                v = v + sumcof(points,self.scatterer{scattererIndex}.getCentre(),kwave,cof(ii),'H');
            end
            
            % insert values into val
            if nargin>5
                val(mask)=v;
            else
                val=v;
            end
            
        end
        
        %-----------------------------------
        % Get useful information about the 
        % matrix.
        %
        % nm(k) - dimension of block for scatterer k
        % nb(k) - dimension of the sub-blocks for each region
        %           for scatterer k
        % aa(k) - indices for scatterer k in the matrix/coefficient
        %           vector
        %-----------------------------------

        function [nm,nb,aa] = matrixDimensions(self)            
            
            % work out dimensions of matrix 
            nm = zeros(1,self.num_scatterers);
            nb = zeros(1,self.num_scatterers);
            for scatk=1:self.num_scatterers
                num_vec = 2 * self.scatterer{scatk}.getNumRegion();
                nb(scatk) = 2*self.nmax(scatk)+1;
                nm(scatk) = num_vec * nb(scatk);                
            end
            
            % setup indices into A for each scatterer
            tmp = cumsum([0,nm]);
            for scatk = 1:self.num_scatterers
                aa{scatk} = tmp(scatk)+1:tmp(scatk+1);
            end
            
        end
        
    end
    
end