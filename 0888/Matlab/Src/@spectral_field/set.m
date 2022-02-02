function f = set(f,varargin)
%SPECTRAL_FIELD set method
%Field on a Gaussian Grid with spherical harmonic coefficients
%  Input:  f= spectral_field([name],[GG],[gp])
%          varargin - ['name',name], ['gp',gp], ['sc', sc]
%  The spectral_field data structure contains
%    f.name - name of the field
%    f.G   -  the gauss_grid on which the field is defined, 
%             with spherical harmonics, etc..
%             Note: G.gg  -  the Gaussian grid (gauss_grid(nj)) ni x nj 
%             oriented north to south (ni = 2*nj)
%    f.gp  -  grid point values of the field (ni x nj) lon-lat
%    f.sc  -  spectral coefficients in a triangular truncations
%             with mtrunc=(2*nj -1 )/3 and kk=mm=nn=mtrunc
% Methods:  get, set, shtrans, shtrana

if (isa(f,'spectral_field'))
     property_argin = varargin;
     while length(property_argin) >= 2,
	prop = property_argin{1};
	val  = property_argin{2};
	property_argin = property_argin(3:end);
	switch prop
	case 'name'
	  f.name = val;
	case 'gp'
	  f.gp = val;
	case 'sc'
	  f.sc = val;
	otherwise
	  error('Set spectral_field method called with unknown property')
        end
     end
else
 error('Call to set spectral_field class method with wrong type ')
end
