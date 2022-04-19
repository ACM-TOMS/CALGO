%% pntsGeometry3D
% Coordinates of points for given geometries in three dimensional space. 
% At the moment only points are arranged on a sphere or have to be defined
% manually.
%
%% Syntax
%
%  [Pnts,Nb,dS] = pntsGeometry3D(seti,type,varargin)
%
%  [Pnts,Nb,dS] = pntsGeometry3D(seti,'manually',Pnts);
%  [Pnts,Nb,dS] = pntsGeometry3D(seti,'manually',Pnts,dispDepth);
%
%  [Pnts,Nb,dS] = pntsGeometry3D(seti,'sphereLatLon',Nb,rad);
%  [Pnts,Nb,dS] = pntsGeometry3D(seti,'sphereLatLon',Nb,rad,dispDepth);
%
%  [Pnts,Nb,dS] = pntsGeometry3D(seti,'sphereFibo',Nb,rad);
%  [Pnts,Nb,dS] = pntsGeometry3D(seti,'sphereFibo',Nb,rad,dispDepth);
%
%
%% Description
%
% |[Pnts,Nb,dS] = pntsGeometry3D(seti,type,varargin)| computes the 
% coordinates in |Pnts| and serves their number |Nb| and an approximation 
% of the infinitesimal element of closed contour in dS. The geometry type 
% is defined in |type| and the details in |varargin| depending on the type.
% |seti| is only used to check that |seti.dim = 3|.
%
% |[Pnts,Nb,dS] = pntsGeometry(seti,'manually',Pnts) counts the number |Nb| 
% of input points |Pnts| and sets |dS| to 1. Input and output |Pnts| is the same.
%
% |[Pnts,Nb,dS] = pntsGeometry3D(seti,'sphereLatLon',Nb,rad)| generates 
% coordinates |Pnts| on a sphere. The number of |Nb| points are arranged on 
% a sphere with radius |rad| with latitude-longitude lattice. 
% Distance |dS| is computed. Note that this in general does not result in a
% equidistant distribution and a fewer number of points is arranged than
% originally given as input. Therefore, we recommend to use the method
% |'sphereFibo'|.
%
% |[Pnts,Nb,dS] = pntsGeometry3D(seti,'sphereFibo',Nb,rad)| generates 
% coordinates |Pnts| on a sphere. The number of |Nb| points are arranged on 
% a sphere with radius |rad| with Fibonacci lattice, see [1,&nbsp;Sec.&nbsp;3.1]. 
% Distance |dS| is computed. This results in a nearly equidistant
% distribution. In case of an even number |Nb|, only |Nb|-1 points are
% arranged; if the input number |Nb| is an odd number nothing will be
% subtracted. Therefore, we recommend to use the method |'sphereFibo'|.
%
% |[Pnts,Nb,dS] = setData(seti,'manually',Pnts,dispDepth)|, 
% |[Pnts,Nb,dS] = pntsGeometry3D(seti,'sphereLatLon',Nb,rad,dispDepth)|,
% and |[Pnts,Nb,dS] = pntsGeometry3D(seti,'sphereFibo',Nb,rad,dispDepth)|
% do the same as without the argument |dispDepth|, but allow to control 
% the depth of displayed messages by |dispDepth|.
%
%% Examples
%
% *Example 1: manually*
%
%   seti.dim = 3;
%   Pnts = [2 0 0 2; 
%           0 2 0 1; 
%           0 0 2 2];
%   [Pnts,Nb,dS] = pntsGeometry3D(seti,'manually',Pnts);
%   figure(101); scatter3(Pnts(1,:),Pnts(2,:),Pnts(3,:),'filled');
%
% <<../extGraph/pntsGeometry3D_fig101.png>>
%
% *Example 2: sphere with latitude-longitude lattice*
%
% Arrange 50 points on a sphere with radius 5. Because latitude-longitude
% lattice actually only 32 points are arranged.
%
%   seti.dim = 3;
%   Nb = 50;
%   rad = 5;
%   [Pnts,Nb,dS] = pntsGeometry3D(seti,'sphereLatLon',Nb,rad);
%   figure(102); [sx,sy,sz] = sphere; r = rad; sx = sx*r; sy = sy*r; sz = sz*r; surf(sx,sy,sz);
%   hold on; scatter3(Pnts(1,:),Pnts(2,:),Pnts(3,:),'filled'); axis square; hold off;
%
% <<../extGraph/pntsGeometry3D_fig102.png>>
%
% *Example 3: sphere with Fibonacci lattice*
%
% Arrange 50 points on a sphere with radius 5. Because 50 is an even
% number, the method with Fibonacci lattice acutally arranges only 49
% points.
%
%   seti.dim = 3;
%   Nb = 50;
%   rad = 5;
%   [Pnts,Nb,dS] = pntsGeometry3D(seti,'sphereFibo',Nb,rad);
%   figure(103); [sx,sy,sz] = sphere; r = rad; sx = sx*r; sy = sy*r; sz = sz*r; surf(sx,sy,sz);
%   hold on; scatter3(Pnts(1,:),Pnts(2,:),Pnts(3,:),'filled'); axis square; hold off;
%
% <<../extGraph/pntsGeometry3D_fig103.png>>
%
%% Input Arguments
%
% * |seti| : structure array which has to contain seti.dim = 3 because three
% dimensional space.
% * |type| : string with type of geometry: 'sphereLatLon', 'sphereFibo'
% (recommended).
% * |varargin| : various input arguments depending on choosen type.
%
% *varargin in case of type = 'manually'*
%
% * |Pnts| : coordinates of points (real matrix of size 3 x Nb).
%
% *varargin in case of type = 'sphereLatLon' and 'sphereFibo'*
%
% * |Nb|          :   Number of points
% * |rad|         :   radius of sphere
%
% *Optional Input Argument*
%
% This argument must be the last if it is used.
%
% * dispDepth   : Depth of displayes messages (0: no, 1 or greater: yes).
%
%% Output Arguments
%
% * |Pnts| : Coordinates of points (real matrix of size 3 x |Nb|).
% * |NbOut| : Number of points (may differ from input |Nb|).
% * |dS| : Approximation of the infinitesimal element of closed contour 
% with control points. A closed contour does not make sense in the cases 
% 'manually' then it is set to 1.
%
%% More About
%
% * The distribution of the points on a sphere with Fibonacci lattice is described in [1,&nbsp;Sec.&nbsp;3.1].
%
%% References
%
% * [1] &Aacute;lvaro Gonz&aacute;lez. Measurement of areas on a sphere using fibonacci and latitude?longitude lattices. _Mathematical Geosciences_, 42(1):49, 2010.
%
%% See Also
%
% * <dS3D.html>
% * <expSetup.html>
%
%% Code: pntsGeometry3D
%

function [Pnts,Nb,dS] = pntsGeometry3D(seti,type,varargin)

if seti.dim ~= 3
    error('pntsGeometry3D needs 3D.')
end

if nargin == 4 && strcmp(type,'manually')
    dispDepth = varargin{2};
else
    dispDepth = 0;
end

%%
% * Split varargin in dependence of parameter type*

if strcmp(type,'manually')
    %   pntsGeometry3D(seti,'manually',Pnts)
    Pnts = varargin{1};
    Nb = size(Pnts,2);
    if nargin == 4
        dispDepth = varargin{2};
    else
        dispDepth = 0;
    end

elseif strcmp(type,'sphereLatLon')
    %   pntsGeometry(seti,'sphereLatLon',Nb,rad)
    Nb = varargin{1};
    rad = varargin{2};
    Pnts = sphereLatLon(Nb,rad);
    if nargin == 5
        dispDepth = varargin{3};
    else
        dispDepth = 0;
    end

elseif strcmp(type,'sphereFibo')
    %   pntsGeometry(seti,'sphereFibo',Nb,rad)
    Nb = varargin{1};
    rad = varargin{2};
    Pnts = sphereFibo(Nb,rad);
    if nargin == 5
        dispDepth = varargin{3};
    else
        dispDepth = 0;
    end
end

% *Compute Nb and dS*

NbOut = size(Pnts,2); % in some cases (especially sphereLatLon, but also 
% in case of an even number sphereFibo) NbOut is not Nb.

if Nb ~= NbOut && dispDepth >= 1
    fprintf('   Input was Nb = %d, but set "Nb" to %d (e.g. because surface mesh of S^2 has %d points).\n',Nb, NbOut,NbOut);
end
Nb = NbOut;

if strcmp(type,'manually')
    dS = 1;
elseif strcmp(type,'sphereLatLon') || strcmp(type,'sphereFibo')
    % normalized surface measure
    dS = dS3D(Nb, rad);
else
    error('pntsGeometry3D.m: no value type to compute dS.');
end

end

%% Code: Subfunctions

%%
% * sphereLatLon*

function Pnts = sphereLatLon(Nb,rad)
%% Arrange points on a sphere with latitude-longitude lattice
%
% Set incident points for a ball in 3D.
%
% Input Arguments:
% * Nb    : number of points (input)
% * rad   : radius of the sphere
%
% Output Arguments:
% * Pnts  : Coordinates of the points.
% * NbOut : number of points (really generated).
%

% incident fields
NbTemp = floor(sqrt(Nb+2));
[pntsX, pntsY, pntsZ] = sphere(NbTemp-1);

% extract multiple points
pntsX = pntsX(:).'; pntsY = pntsY(:).'; pntsZ = pntsZ(:).';
Pnts = [pntsX(1); pntsY(1); pntsZ(1)];
for jj = 2:length(pntsX)
    [~,aux] = size(Pnts);
    dummy = 1;
    for kk = 1:aux
        if norm(Pnts(:,kk)-[pntsX(jj); pntsY(jj); pntsZ(jj)]) < 1e-10
            dummy = 0;
        end
    end
    if dummy == 1
        Pnts(:,aux+1) = [pntsX(jj); pntsY(jj); pntsZ(jj)];
    end
end

[~,NbOut] = size(Pnts); % no output of NbOut, but if below needs it.

Pnts = rad*Pnts;

if NbOut == 1
    % workaround: code above does not work for Nb == 1.
    if dispDepth >= 1
        disp('pntsGeometry3D.m : As Nb == 1 => Set Pnts = [1; 0; 0]');
    end
    Pnts = [1; 0; 0];
    % NbOut = 1; % no output here
    % dS = 1; % no output here
end

end

%%
% *sphereFibo*

function Pnts = sphereFibo(Nb,rad)
% Arrange points on a sphere with Fibonacci lattice
%
% * If Nb is odd, then really Nb points are generated.
% * If Nb is even, Nb-1 points are generated.

N = floor((Nb-1)/2); 

%%
% *Geographical coordinates in degrees*
%
% The following follows the pseudo code in Sec. 3.1 in:
% &Aacute;lvaro Gonz&aacute;lez. Measurement of areas on a sphere using fibonacci and latitude?longitude lattices. _Mathematical Geosciences_, 42(1):49, 2009.
% 

NN = 2*N+1; % resulting number of points...
lat = zeros(1,NN); % latitude (North/South)
lon = zeros(1,NN); % longitude (West/East)
phi = (1+sqrt(5))/2; % phi: golden ratio
for i = -N:N
    j = i+N+1; % shifted i because MATLAB can not deal with negative index (or 0)
    lat(j) = asin(2*i/(2*N+1))*180/pi; % asin is arcsin
    % asin is arcsin, but of course needs [-1,1] input for a real output.
    lon(j) = mod(i,phi)*360/phi;
    if lon(j) < -180
        lon(j) = 360 + lon(j);
    elseif lon(j) > 180
        lon(j) = lon(j) - 360;
    end
end

%%
% *Coordinates in cartesian system*
%
% Transforms geographical coordinates (they are spherical coordinates) into
% cartesian. Note that |lat| and |lon| are in degrees and not in rad.
%
% Transform them to rad
latrad = deg2rad(lat);
lonrad = deg2rad(lon);
%
% Transform them into cartesian system
r = 1;
[x,y,z] = sph2cart(lonrad,latrad,r);
% n = size(x,2); % number of points: no output here

Pnts = zeros(1,NN);
Pnts(1,:) = x; Pnts(2,:) = y; Pnts(3,:) = z;

Pnts = rad*Pnts;

end
