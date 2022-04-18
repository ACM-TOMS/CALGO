function g = test(f)
% TEST: Spectral synthesis and analysis test program for the spectral_field class
% Usage:  test(f)
% Input:   f is a spectral_field instantiated before call to test().
% Output:  on return there is no change to the object f
%          plots are produced of test field, energy spectrum, inverted field,
%          and the difference between untransformed and transformed test field.
%--------------------------------------------------------
%   <Description>
% Field on a Gaussian Grid with spherical harmonic coefficients
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
% Methods:  get, set, shtrans, shtrana, test, del2, del2inv
%--------------------------------------------------------
% Example usage:
%  G = gauss_grid('T10',16);
%  f = spectral_field('test field',G);
%  test(f);
%--------------------------------------------------------
% Local variables
%    x - complex spectral coeffecients 
%    nj = number of Gauss latitudes
%    mm, nn,kk are the truncation parameters
%    yg - Gauss points in (-1.1) indexed by j
%    wg - Gauss weights
%    P - associated Legendre functions ordered (j,n,m)
%    f - grid point values
%--------------------------------------------------------
ni = get(f.G,'ni');
nj = get(f.G,'nj');
mm = get(f.G,'mm');
nn = get(f.G,'nn');
kk = get(f.G,'kk');
wg = get(f.G,'wg');
xg = get(f.G,'xg');
yg = get(f.G,'yg');
P  = get(f.G,'P');
%--------------------------------------------------------
%Fill a field array with a test function
%--------------------------------------------------------
gp= zeros(ni,nj);
m=3;n=4;
sgn = (-1)^(n-m);
for i=1:ni
%    emil = cos(m*xg(i));
%    emil = sin(m*xg(i));
    emil = cos(m*xg(i))+ sin(m*xg(i)); 
%    emil = 1.0;
    for j=1:nj/2
        gp(i,j)=emil*P(j,n+1,m+1);
        gp(i,nj+1-j)=sgn*emil*P(j,n+1,m+1);
%        gp(i,j)=emil;
%        gp(i,nj+1-j)=emil;
        
    end
end
f=set(f,'gp',gp)
plot(f, 'Test function plot of Y_4^3')
pause
ftest = f;
%--------------------------------------------------------
%  Compute the forward analysis using a spherical harmonic transform
%--------------------------------------------------------
%profile on;
timanalysis=cputime;
f=shtrana(f);
timanalysis=cputime - timanalysis
%profile report
%--------------------------------------------------------
%  Get the spectral coefficients of the transformed field
%--------------------------------------------------------
x = get(f,'sc');
%x                                 !print out spectral coefficients
%--------------------------------------------------------
%Compute power spectrum
%--------------------------------------------------------
e = zeros(1,nn);
for n=1:nn
    engy=0.0;
    for m=0:n
        engy = engy + real( x(n+1,m+1)*conj(x(n+1,m+1)) );
    end
    e(n) = engy;
end
%--------------------------------------------------------
% plot the energy spectrum
%--------------------------------------------------------
plot(e)
title('Test Energy Spectrum')
pause;
%--------------------------------------------------------
% compute L2 norm of the function
%--------------------------------------------------------
s = l2norm(f)
%--------------------------------------------------------
% Transform back to grid point space for error check
%--------------------------------------------------------
timsynthesis=cputime;
f=shtrans(f);
timsynthesis=cputime - timsynthesis
%--------------------------------------------------------
% Get the grid point values
%--------------------------------------------------------
gp2 = get(f,'gp');
%--------------------------------------------------------
% Get the xg and yg from the gauss_grid and plot
%--------------------------------------------------------
%h = surf(yg,xg,gp2);
plot(f,'Test (shtrana/shtrans) forward/backward transformed field')
pause;
herr = surf(yg,xg,gp-gp2);
title('Test error in forward/backward transformed field')
pause
%--------------------------------------------------------
% Compute the Laplacian(del2) of the test field
%--------------------------------------------------------
g = del2(f);
%gp2 = get(g,'gp');
%hlapl = surf(yg,xg,gp2);
plot(g,'Test (del2) Laplacian of the field ')
pause
%--------------------------------------------------------
% Compute the inverse Laplacian(del2iinv) of the test field
%--------------------------------------------------------
g = del2inv(g);
gp2 = get(g,'gp');
%hlapl = surf(yg,xg,gp2);
plot(g,'Test (del2inv) inverse Laplacian ')
pause
hlapl = surf(yg,xg,gp-gp2);
title('Error (f-del2inv(del2(f)))')
pause
%--------------------------------------------------------
% Compute the spherical Helmholtz equation solution of the test field
%--------------------------------------------------------
c2 = 1.0;
g = helmholtz(c2,g);
gp2 = get(g,'gp');
plot(g, 'Test (helmholtz) Helmholtz solve ')
pause
%--------------------------------------------------------
% Compute the spherical gradient of the test field
%--------------------------------------------------------
g = ftest;
[g1,g2] = grad(g);
plot(g1,'Test (grad1) gradient ')
g1
pause
plot(g2,'Test (grad2) gradient ')
g2
pause
%--------------------------------------------------------
% Compute the spherical (horizontal) curl of the test field
%--------------------------------------------------------
g = curl(g1,g2)      %curl (grad f) == 0 
plot(g,'Test (curl) ')
pause
%--------------------------------------------------------
% Compute the spherical divergence of the test field
%--------------------------------------------------------
g = div(g1,g2); 
%gp= get(g,'gp');
%hcurl = surf(yg,xg,gp);
plot(g,'Test (div) Divergence ')
pause
%--------------------------------------------------------
% Compute the potential vorticity inversion
%--------------------------------------------------------
[U,V] = UVinv(g1,g2);
gpu= get(U,'gp');
gpv= get(V,'gp');
%hcurl = surf(yg,xg,gpu);
plot(U,'Test (UVinv: U) Potential Vorticity Inversion ')
pause
%hcurl = surf(yg,xg,gpv);
plot(V,'Test (UVinv: V) Potential Vorticity Inversion ')
pause

