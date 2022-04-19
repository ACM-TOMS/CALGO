
function A = slt_grid(varargin)
%SLT Creates an extended Grid from an input Grid
%  Usage:  Adv = slt(xg,yg)
%    Once the extended grid is established, with index mappings for
%    halo updates, then particle tracking and interpolation at departure
%    points can be done from input fields.
% Input:  G  - grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  The grid data structure contains
%    ni, nj           the number of internal lons and lats
%    nie, nje         the number of total lons and lats in extend slt_grid
%    xg,yg   -  the Gaussian grid ni x nj oriented north to south (ni = 2*nj)
%  The extended grid for advection contains
%    nix, njx
%    xge,yge  -   the extended grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch nargin
case 0 
% no input arguments so create a default component
 A.xge     = [];
 A.yge     = [];
 A.X       = [];
 A.Y       = [];
 A.Xe       = [];
 A.Ye       = [];
 A.nie    = 0;
 A.nje	  = 0;
 A.ni     = 0;
 A.nj	  = 0;
 A.nxpt   = 0;
 A.nypt   = 0;
 A = class(A,'slt_grid');
case 1
% one input argument so check if it is a component already
 if (isa(varargin{1},'slt_grid'))
   A= varargin{1};
 end
case 2
 xg   = varargin{1};
 yg   = varargin{2};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% form the slt extended grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  -------------------------------------
%  |                                   |
%  |       C                D          |
%  |  +++++++++++++++N++++++++++++++   +1.0
%  |  +                            +   |
%  |  +                            +   |
%  |  +            interior        + periodic
%  |  +                            +   |
%  |  +                            +   |
%  |  +++++++++++++++S++++++++++++++   -1.0 (= sin lat)
%  | 0.0   A                B      2*pi
%  |                                   |
%  |                                   |
%  -------------------------------------
%     +nxpt+1                      +nxpt+nx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 nx = length(xg);
 ny = length(yg);
% Periodic Longitude
 nxpt = 3;  % number of extra points in periodic x extension
 exl = [xg(nx-nxpt+1) xg(nx-nxpt+2) xg(nx-nxpt+3) ]- 2*pi;
 exr = [xg(1) xg(2) xg(3) ]+ 2*pi;
 xge   = [exl xg exr];
% Latitude as \sin \phi in [-1,1]
 nypt = 3;  % number of extra points in pole y extension (skipping the pole)
 eyb = [(-2 - yg(3)) (-2 - yg(2)) (-2 - yg(1)) ];       % bottom (south pole)
 eyt = [(2 - yg(ny)) (2 - yg(ny-1)) (2 - yg(ny-2)) ];   % top    (north pole)
 yge   = [eyb yg eyt];
%Note:  ESMF_COORD_SYSTEM_SPHERICAL  (or lat,lon)
%Note:  ESMF_COORD_ORDER_UNKNOWN  (or lat,lon)
%Note:  ESMF_GRID_TYPE_LATLON
%Note:  ESMF_GRID_TYPE_LATLON
%Fill the coordinate arrays
 A.xge = xge;
 A.yge = yge;
 [A.X A.Y]   = meshgrid(xg,yg);
 [A.Xe A.Ye] = meshgrid(xge,yge);
 A.nie       = length(A.xge);
 A.nje       = length(A.yge);
 A.ni       = length(xg);
 A.nj       = length(yg);
 A.nxpt      = nxpt;
 A.nypt      = nypt;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end slt_fill
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 A = class(A,'slt_grid');
otherwise
 error('Wrong number of input arguments to slt_grid class constructor')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end slt_grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

