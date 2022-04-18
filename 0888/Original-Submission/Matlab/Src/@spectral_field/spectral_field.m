function f = spectral_field(varargin)
%SPECTRAL_FIELD Constructor method
%Field on a Gaussian Grid with spherical harmonic coefficients
%    spectral_field([name],[GG],[gp])
%    The spectral coefficients may be obtained by a spherical
%    harmonic transform of the grid point values.  The transform and the inverse
%    transform are public methods of this class.
%  The data structure contains
%    f.name - (string) name of the field, eg. 'temperature'
%    f.G   -  the gauss_grid on which the field is defined, 
%             with spherical harmonics, etc..
%             Note: G.gg  -  the Gaussian grid (gauss_grid(nj)) ni x nj 
%             oriented north to south (ni = 2*nj)
%    f.gp  -  grid point values of the field (ni x nj) lon-lat
%    f.sc  -  spectral coefficients in a triangular truncations
%             with mtrunc=(2*nj -1 )/3 and kk=mm=nn=mtrunc
% Methods:  get, set, shtrans, shtrana, test, plot

switch nargin
case 0
% no input arguments so create a default component
 f.name  ='';
 f.G     = {};
 f.gp    = [];
 f.sc    = [];
 f       = class(f,'spectral_field');
case 1
% one input argument so check if it is a component already
 if (isa(varargin{1},'spectral_field'))
   f= varargin{1};
 else
 f.name  = varargin{1};
 f.G     = {};
 f.gp    = [];
 f.sc    = [];
 f       = class(f,'spectral_field');
 end
case 2
% two input arguments so initialize name and gauss_grid
 f.name  = varargin{1};
 f.G     = varargin{2}; % gauss_grid type
 ni      = get(f.G,'ni');
 nj      = get(f.G,'nj');
 f.gp    = zeros(ni,nj);
 nn      = get(f.G,'nn');
 mm      = get(f.G,'mm');
 f.sc    = complex(zeros(mm,nn),zeros(mm,nn));
%Note:  ESMF_COORD_SYSTEM_SPHERICAL  (or lat,lon)
%Note:  ESMF_GRID_TYPE_LATLON
 f       = class(f,'spectral_field');
case 3
% two input arguments so initialize name and other variables
 f.name  = varargin{1};
 f.G     = varargin{2} % gauss_grid type
 f.gp    = varargin{3}; % grid point values should be (G.ni,G.nj)
 nn      = get(f.G, 'nn');
 mm      = get(f.G, 'mm');
 f.sc    = complex(zeros(mm,nn),zeros(mm,nn));
%Note:  ESMF_COORD_SYSTEM_SPHERICAL  (or lat,lon)
%Note:  ESMF_GRID_TYPE_LATLON
 f       = class(f,'spectral_field');
otherwise
 error('Wrong number of input arguments to spectral_field class constructor')
end
