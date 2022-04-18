function G = gauss_grid(varargin)
%Gauss_grid Creates a Gaussian Grid with Gauss weights and the
%    associated Legendre functions evaluated on this grid.
% Input:  name, nj
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  The data structure contains
%    ni, nj           the number of lons and lats
%    mm,nn,kk         with mtrunc=(2*nj -1 )/3 and kk=mm=nn=mtrunc, settings for
%                     a triangular spectral truncation
%    gg   -  the Gaussian grid ni x nj oriented north to south (ni = 2*nj)
%    wg   -  Gauss weights, nj
%    P    -  associated Legendre functions evaluated at the latitudes
%    H    -  derivative of associated Legendre functions evaluated at the latitudes
%  Note on usage:  this is primarily used for computations of spectral fields.
%                  The spectral transform of a field over this grid must pass the
%                  handle of this object to the methods.

switch nargin
case 0 
% no input arguments so create a default component
 G.name  ='';
 G.ni     = 0;
 G.nj	  = 0;
 G.mm     = 0;
 G.nn     = 0;
 G.kk     = 0;
 G.radius = 1.0;   % default unit sphere
 G.xg     = [];
 G.yg     = [];
 G.wg     = [];
 G.P      = [];
 G.H      = [];
% G.gg     = [];
 G = class(G,'gauss_grid');
case 1
% one input argument so check if it is a component already
 if (isa(varargin{1},'gauss_grid'))
   G= varargin{1};
 else
   G.name  = varargin{1};
   G.ni     = 0;
   G.nj	    = 0;
   G.mm     = 0;
   G.nn     = 0;
   G.kk     = 0;
   G.radius = 1.0;   % default unit sphere
   G.xg     = [];
   G.yg     = [];
   G.wg     = [];
   G.P      = [];
   G.H      = [];
%   G.gg     = [];
   G = class(G,'gauss_grid');
 end
case 2
% two input arguments so initialize name and other variables
 G.name  = varargin{1};
 G.nj = varargin{2}; % nj is number of latitudes
 G.ni=2*G.nj;               % number of longitudes
 mtrunc = floor((2*G.nj -1)/3);
 G.mm=mtrunc; G.nn=mtrunc; G.kk=mtrunc; % triangular spectral truncation (nn=mm=kk)
 G.radius = 1.0;   % default unit sphere (can be set later on)
 [G.xg,G.yg,G.wg,G.P,G.H] = shtraninit(G.ni,G.nj,G.mm,G.nn,G.kk);
%Note:  ESMF_COORD_SYSTEM_SPHERICAL  (or lat,lon)
%Note:  ESMF_GRID_TYPE_LATLON
%Fill the coordinate arrays
% G.gg = meshgrid(G.xg,G.yg);
 G = class(G,'gauss_grid');
otherwise
 error('Wrong number of input arguments to gauss_grid class constructor')
end

