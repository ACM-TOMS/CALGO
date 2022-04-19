more off
echo on
clc
% This script demonstrates Faber polynomials.

% Faber polynomials are defined via a conformal map f from the
% simply connected exterior of a bounded region to the exterior
% of the unit disk, fixing the point at infinity.  The nth Faber 
% polynomial is the polynomial part of the Laurent expansion of f^n
% at infinity.  The Faber polynomials reduce to Chebyshev polynomials
% when the region is an interval and Taylor polynomials when it is a
% disk.  If the region is bounded by a polygon, the Schwarz-Christoffel
% exterior map can be used to compute Faber polynomial coefficients.

% If you are using MATLAB for MS Windows, be sure you have selected the
% "Enable background process" item on the "Options" menu before proceeding.

pause    % Strike any key to begin (Ctrl-C to abort)
figure(gcf)
hold off
clc
% Use the mouse to draw a polygon.  Be sure to put vertices in
% clockwise order, and use only finite vertices.

[w,beta] = drawpoly; 

hold on
axis(axis)

pause    % Strike any key to compute Faber polynomial coefficients

[z,c] = deparam(w,beta);

F = faber(20,w,beta,z,c);

pause    % Strike any key to continue


% Because the Faber polynomials approximate a function having unit
% modulus on the polygon, the lemniscates {z: |p(z)|=1} for Faber
% polynomials p will approximate the polygon.

lim = axis;
[X,Y] = meshgrid(linspace(lim(1),lim(2),40),linspace(lim(3),lim(4),40));
h = line(NaN,NaN);
for m = 4:4:16
  delete(h)
  Z = abs(polyval(F(1:m+1,m+1),X+i*Y));          
  [con,h] = contour(X,Y,Z,[1,NaN],'c');
  title(['degree of Faber polynomial = ',int2str(m)])
  disp(' ')
  disp([blanks(m/4),'   Strike any key to continue'])
  pause
end


% Another way to see this is to look at |p(z)| for z on the polygon.
% Here we choose Fejer points (images of roots of unity) for z.

zp = demap(exp(i*linspace(0,2*pi,100)),w,beta,z,c);

hold off
for m = 4:4:16
  plot(1:100,abs(polyval(F(1:m+1,m+1),zp)));          
  title(['degree of Faber polynomial = ',int2str(m)])
  disp(' ')
  disp([blanks(m/4),'   Strike any key to continue'])
  pause
end

echo off    % End of demo
title(' ')
hold off


