function val = get(f,prop)
%SPECTRAL_FIELD get method
%Field on a Gaussian Grid with spherical harmonic coefficients
%  Input:  f= spectral_field([name],[GG],[gp])
%          prop - ['name'], ['gp'], ['sc']
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
	switch prop
	case 'name'
	  val = f.name;
	case 'gp'
	  val = f.gp;
	case 'sc'
	  val = f.sc;
	otherwise
	  error('Get spectral_field method called with unknown property')
        end
else
 error('Call to get spectral_field class method with wrong type ')
end
