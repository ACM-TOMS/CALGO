function val = get(G,prop_name)
%GET method for Gauss_grid class 
% Usage:  val = get(G,prop_name);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  The gauss_grid data structure contains
%    ni, nj           the number of lons and lats
%    mm,nn,kk         with mtrunc=(2*nj -1 )/3 and kk=mm=nn=mtrunc, settings for
%                     a triangular spectral truncation
%    xg, yg  - the array of lons and lats
%    radius  - of the sphere
%    wg   -  Gauss weights, nj
%    P    -  associated Legendre functions evaluated at the latitudes
%    H    -  derivative of associated Legendre functions evaluated at the latitudes
%    gg   -  the Gaussian grid ni x nj oriented north to south (ni = 2*nj)
%  Note on usage:  this class is primarily used for computations of spectral fields.
%                  The spectral transform of a field over this grid must pass the
%                  handle of this object to the methods.
% no input arguments so create a default component
 switch prop_name
 case 'nj'
	val = G.nj;
 case 'ni'
	val = double(G.ni);
 case 'mm'
	val = G.mm;
 case 'nn'
	val = G.nn;
 case 'kk'
	val = G.kk;
 case 'radius'
	val = G.radius;
 case 'xg'
	val = G.xg;
 case 'yg'
	val = G.yg;
 case 'wg'
	val = G.wg;
 case 'gg'
%	val = G.gg;
 case 'P'
	val = G.P;
 case 'H'
	val = G.H;
otherwise
	error([prop_name , 'Is not a gauss_grid property'])
end

