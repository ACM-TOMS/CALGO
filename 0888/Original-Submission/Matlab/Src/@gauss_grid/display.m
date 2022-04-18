function display(G)
%Display method for Gauss_grid class 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  The gauss_grid data structure contains
%    ni, nj           the number of lons and lats
%    mm,nn,kk         with mtrunc=(2*nj -1 )/3 and kk=mm=nn=mtrunc, settings for
%                     a triangular spectral truncation
%    xg, yg  - the array of lons and lats
%    wg   -  Gauss weights, nj
%    P    -  associated Legendre functions evaluated at the latitudes
%    H    -  derivative of associated Legendre functions evaluated at the latitudes
%    gg   -  the Gaussian grid ni x nj oriented north to south (ni = 2*nj)
%  Note on usage:  this class is primarily used for computations of spectral fields.
%                  The spectral transform of a field over this grid must pass the
%                  handle of this object to the methods.
% no input arguments so create a default component
 disp(' ');
 disp([inputname(1) ' '])
 stg=sprintf('name: %s',G.name);
 disp(stg);
 stg=sprintf('nj: %5i',G.nj);
 disp(stg);
 stg=sprintf('ni: %5i',G.ni);
 disp(stg);
 stg=sprintf('mm: %5i',G.mm);
 disp(stg);
 stg=sprintf('nn: %5i',G.nn);
 disp(stg);
 stg=sprintf('kk: %5i',G.kk);
 disp(stg);
 stg=sprintf('radius: %8g',G.radius);
 disp(stg);

