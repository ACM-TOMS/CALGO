function G = set(G,prop_name,val)
%SET method for Gauss_grid class 
% Usage:  G = set(G,prop_name,val);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  The gauss_grid data structure contains
%    ni, nj           the number of lons and lats
%    mm,nn,kk         with mtrunc=(2*nj -1 )/3 and kk=mm=nn=mtrunc, settings for
%                     a triangular spectral truncation
%    radius           of the sphere
%    xg, yg  - the array of lons and lats
%    wg   -  Gauss weights, nj
%    P    -  associated Legendre functions evaluated at the latitudes
%    H    -  derivative of associated Legendre functions evaluated at the latitudes
%    gg   -  the Gaussian grid ni x nj oriented north to south (ni = 2*nj)
%  Note on usage:  this class is primarily used for computations of spectral fields.
%                  The spectral transform of a field over this grid must pass the
%                  handle of this object to the methods.
% no input arguments so create a default component
 switch prop_name
 case 'radius'
	G.radius= val;
otherwise
	error([prop_name , 'is not set-able'])
end

