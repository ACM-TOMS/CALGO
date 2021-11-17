
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

classdef mfsHexagonSolver < mfsSolver
    
    properties
        cornern
    end
    
    methods
        
        %-----------------------------------------
        % constructor
        %-----------------------------------------
        
        function self = mfsHexagonSolver(kwave,incidentField,n,tau,m)
            
            % set defaults for MFS parameters
            if nargin < 6
                tau = 5e-2;
            end
            
            if nargin < 7
                m = 2*n;
            end
            
            %  call parent constructor
            self = self@mfsSolver(kwave,incidentField,n,tau,m);            

            % set number of corner basis functions
            self.cornern = n;
            
        end
            
        %-----------------------------------------
        % setup
        %-----------------------------------------
        
        function setup(self)

            % setup structure with the parameters for the MFS method
            opts= struct('eta',self.kwave,'fast',2,'nmultiplier',2,'tau',self.tau);
            
            % setup structure with parameters for the corner basis
            % functions
            nuopts=struct('type','s','cornermultipliers',[0 0 1 0 0],'rescale_rad',1);

            % set number of vertices
            V = 6;
            
            % set radius of polygon
            r0 = 1;
            
            % set radius of bounding circle
            r1 = 2;

            % compute polar angles of the radial edges of the artifical
            % domains
            t = 2*pi*((0:V-1)-0.5)/V;
            
            % compute corners of artificial domains that lie on the
            % bounding circle
            r = r1 * exp(1i*t);
            
            % compute the corners of the polygon 
            a = r0 * exp(1i*2*pi*(0:V)/V);

            % compute the midpoints of the faces
            c = 0.5*(a(2:end)+a(1:end-1));
            
            % set the radius for the MFS basis
            R = 0.8*r1;
            
            % create straight line segments that are used to construct the 
            % first artifical domain 
            s = segment.polyseglist(self.m,[r(2) c(1) a(1) c(V) r(1)]);
            
            % add the boundary that lies on the bounding circle
            s = [s(1:3) segment(3*self.m,[0 r1 t(1) t(2)])];

            % create a temporary copy that we can rotate to make the other
            % artificial domains
            tmp = [s(1:3) segment(3*self.m,[0 r1 t(1) t(2)])];
            
            % rotate the temporary copy to make the other artificial
            % domains
            for k=1:V-1
                s = [s rotate(tmp,2*pi*k/V)];
            end
            
            % collect segments that make up all the  artificial boundaries
            sart = s(interlace(1+4*(0:V-1),4+4*(0:V-1)));
            
            % collect segments that make up the artificial outer boundary
            sext = s([4+4*(V-1:-1:0)]);
            
            % collect segments that make the polygon
            polygon = s(interlace(2+4*(V-1:-1:0),3+4*(V-1:-1:0)));

            % setup domains for the artifical domains
            for k=1:V
                ii = [1 2 3 -3 4] + 4*(k-1);
                ii = mod(ii-1,4*V)+1;
                d(k)=domain(s(ii),[1 1 1 -1 1]);
            end
            
            % setup a domain for the polygon
            dpoly = domain([],[],polygon,1);
            
            % setup the exterior domain
            ext = domain([], [], sext, -1);
            
            % set transmission BC on the artificial boundaries
            sart.setmatch([-self.kwave self.kwave],[1 -1]);

            % add the corner basis functions
            for k=1:V
                d(k).addcornerbases(self.cornern,nuopts);
            end
            
            % setup MFS basis for the exterior
            Z=@(t) R*exp(2i*pi*t);
            Zp=@(t) 2i*pi*R*exp(2i*pi*t);            
            ext.addmfsbasis({Z,Zp},self.n,opts);
            
            % initialise the scattering problem
            self.scatteringObject = scattering(ext,d);

        end
                
    end % end methods
    
end