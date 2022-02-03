
function s = l2norm(f) 
% L2NORM - spherical L2 norm of the spectral function
% Usage: s = l2norm(f)
% Compute the L2 norm of the function from the grid point values using the spectral
% transform and summing the spectral coefficients (Parseval).  The return
% value is a scalar.
%--------------------------------------------------------
%  Input: 
%    f - spectral_field class object
%        In particular:
%           f.G  - field gauss_grid
%           f.gp - field grid point values
%           f.sc - field spectral coefficients
%  Output:
%    g - spectral_field object = del2(f)
%  Local
%    xf - matrix of Fourier coefficients ordered (m,j)
%--------------------------------------------------------
% Written by John Drake
% Method of spectral_field class: Sept 2005
%--------------------------------------------------------
mm = get(f.G,'mm');
nn = get(f.G,'nn');
kk = get(f.G,'kk');
aradius= get(f.G,'radius');
%--------------------------------------------------------
g = shtrana(f);    % transform the field to spectral space
sc = get(g,'sc');   % get the spectral coeffients
% modify the spectral coefficients
 s = 0.0;
for n = 1:nn
    for m = 0:n
	  s = s + real( sc(n+1,m+1)*conj(sc(n+1,m+1)) );  %is correct
    end
end
s = sqrt(s);
% all done  - f is unchanged but may have computed spectral coefficients.
