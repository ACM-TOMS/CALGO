
function g = helmholtz(alpha, f) 
% HELMHOLTZ - Solve the Helmholtz equation for a spectral_grid field.
% Usage: g = helmholtz(alpha,f)
%     equation:  (Ig + alpha Lapl(g)) = f  
%     alpha (real or complex, possibly negative)  .
% Compute the Helmholtz equation solution 
%  from the grid point values of the right hand side f using the spectral
% transform and modifying the spectral coefficients.  The grid point
% values of the solution g are returned from a spherical harmonic transform synthesis
%--------------------------------------------------------
%  Input: 
%    f - spectral_field class object
%        In particular:
%           f.G  - field gauss_grid
%           f.gp - field grid point values
%           f.sc - field spectral coefficients
%  Output:
%    g - spectral_field object = helmoltz solution with rhs (f)
%  Local
%    sc - matrix of spectral coefficients ordered (n,m)
%--------------------------------------------------------
% Written by John Drake
% Method of spectral_field class: May 2006
%--------------------------------------------------------
mm = get(f.G,'mm');
nn = get(f.G,'nn');
kk = get(f.G,'kk');
aradius = get(f.G,'radius');
%--------------------------------------------------------
g = shtrana(f);    % transform the field to spectral space
sc = get(g,'sc');  % get the spectral coeffients
% modify the spectral coefficients for Helmholtz solution
%
for n = 1:nn
   	m = 0:n;
	sc(n+1,m+1) = sc(n+1,m+1)/(1.0 - alpha*n*(n+1)/aradius^2) ;
end
%  Truncate the approximation for aliasing
sc(nn+1,:) = complex(0.0,0.0);
g = set(g,'sc',sc);% set the new coefficients
g = shtrans(g);    % inverse transform to get new grid point values 
% all done  - f is unchanged, but new spectral_field g created.
