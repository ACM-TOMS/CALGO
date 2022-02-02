function display(f)
%SPECTRAL_FIELD display method
%Field on a Gaussian Grid with spherical harmonic coefficients
%  The data structure contains
%    f.name - name of field
%    f.G   -  the gauss_grid on which the field is defined, 
%             with spherical harmonics, etc..
%             Note: G.gg  -  the Gaussian grid (gauss_grid(nj)) ni x nj 
%             oriented north to south (ni = 2*nj)
%    f.gp  -  grid point values of the field (ni x nj) lon-lat
%    f.sc  -  spectral coefficients in a triangular truncations
%             with mtrunc=(2*nj -1 )/3 and kk=mm=nn=mtrunc
 f= shtrana(f);   % compute the spectral coefficients to display
 disp(' ');
 disp([inputname(1) ' '])
 stg=sprintf('spectral_field name: %s',f.name);
 disp(stg);
 stg=sprintf('gauss_grid: ');
 disp(stg);
 display(f.G);
 stg=sprintf('Field grid point values, gp: ');   % real numbers
 disp(stg);
 f.gp
 stg=sprintf('Spectral coefficients, sc: ');  % complex numbers
 disp(stg);
 f.sc

