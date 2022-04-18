more off
echo on
clc
% This script demonstrates mapping to elongated polygons with the
% Schwarz-Christoffel Toolbox.

pause     % Strike any key to begin (Ctrl-C to abort)


% We begin with a demonstration of "crowding."

pause     % Strike any key to continue

% Here is a moderately elongated region.
w = [3+2.4i; .6-.4i; -1-.4i; -3+2i; -3-2i; -1-.8i; .6-.8i; 3-2i];
beta = scangle(w);

figure(gcf)
hold off
plotpoly(w,beta)

pause     % Strike any key to continue

% Solve the parameter problem for the half-plane.
[x,c] = hpparam(w,beta);

pause     % Strike any key to display results

hpdisp(w,beta,x,c)

pause     % Strike any key to continue

% Notice how close together four of the prevertices are.  
% Although we estimate an aspect ratio of 6/.4 = 15, the 
% prevertices differ by about 1e-7, or roughly exp(-15).
% The exponential proximity of the prevertices is commonly
% known as the crowding phenomenon.  The phenomenon occurs
% for the disk as well as the half-plane (and for exterior
% maps, when the exterior region is elongated).

pause     % Strike any key to continue
clc
% If the polygon above had been much more elongated, the 
% prevertices would have been indistinguishable in double
% precision.  Even before this point, the solution of the
% parameter problem can become extremely difficult.

% While sometimes an elongated polygon can be subdivided 
% and then mapped, a more elegant solution is to use a more
% appropriate fundamental domain.  In the important case of 
% a region which is elongated in only one direction, a natural
% choice is a rectangle.  

pause     % Strike any key to continue
clc
% Let's map the same polygon to a rectangle.  The corners of 
% the rectangle should map to the outermost vertices of the
% polygon. 

[z,c,L]=rparam(w,beta,[1,4,5,8]);

rdisp(w,beta,z,c,L)

pause     % Strike any key to continue

% The corners of the rectangle are found to be about +-pi/2
% and +-pi/2 + 23i. The conformal modulus of the polygon,
% which is the aspect ratio of the rectangle, is determined 
% (as part of the solution) to be about 7.3. 

pause     % Strike any key to continue

% Here's a plot of the images of 6 vertical and 12 horizontal
% lines.

rplot(w,beta,z,c,L,6,12)

pause     % Strike any key to continue
clc
% Here's another example of a rectangle map.  The solution is
% given, just to save time. 
w = [-3 + 1.5i;-2 + 1.5i;1.5 + 1.5i;1.5 + 0.5i;-2 + 0.5i;-2 - 2i;
   3 - 2i;3 - 1i;-1.5 - 1i;-1.5;2;2 + 2.5i;-3 + 2.5i];
 
beta = scangle(w);

z = [
 -0.75095889852766                    
  1.57079632679490                    
  1.57079632679490+13.62792581256858i
  1.57079632679490+21.47908439726175i
  1.57079632679490+42.08393917979348i
  1.57079632679490+51.32236373477432i
  1.57079632679490+65.40343136968455i
 -1.57079632679490+65.40343136968455i
 -1.57079632679490+49.93558255696847i
 -1.57079632679490+42.08442599909570i
 -1.57079632679490+21.47957121760661i
 -1.57079632679490+12.24115478564419i
 -1.57079632679490 ];

c = 1.549060450542251e+13 + 1.549060450542251e+13i;

L = 20.37728759500835;
 
% In the plot, notice how one rectangle corner is mapped to a
% trivial vertex---one located in the middle of a side.

pause     % Strike any key for plot

rplot(w,beta,z,c,L,6,12)
 
pause     % Strike any key to continue
clc
% Another choice for the fundamental domain for an elongated
% polygon is the strip 0 <= Im z <= 1.  This is especially
% appropriate when the target region is a polygonal channel,
% such as you might encounter in a fluids problem.

pause     % Strike any key to continue

% Here's a simple example:

w = [-2-i; -2-2i; -2i; 2-i; Inf; 2.5; Inf];
beta = [.5; -.5; -atan(1/2)/pi; atan(1/2)/pi; -1.2; .2; -1];
plotpoly(w,beta)

pause     % Strike any key to solve the parameter problem

[z,c] = stparam(w,beta,[7,5]);

pause     % Strike any key to see results

stdisp(w,beta,z,c)

stplot(w,beta,z,c)

pause     % Strike any key to continue
clc
% We close with another strip map.  This time, one end of the
% strip will map to a finite vertex.  The result is as if a
% sink or source were placed at this vertex.

w = [3-.5i; 2+1.5i; .5+.5i; -1+.5i; Inf; -1; .5; 2-1.5i];
beta = scangle(w);
beta(4:6) = [.2;-1.4;.2];

[z,c] = stparam(w,beta,[5,1]);
stplot(w,beta,z,c,12,8)

echo off 	% End of demo

