more off
echo on
clc
% This script demonstrates the basic capabilities of the
% Schwarz-Christoffel Toolbox.

pause     % Strike any key to begin (Ctrl-C to abort)


% We begin with an L-shaped region.

pause     % Strike any key to continue

w = [i; -1+i; -1-i; 1-i; 1; 0];
beta = scangle(w);

figure(gcf)
hold off
plotpoly(w,beta)

pause     % Strike any key to continue

% Now we solve the parameter problem for the half-plane.
[x,c] = hpparam(w,beta);

pause     % Strike any key to display results

hpdisp(w,beta,x,c)

pause     % Strike any key to continue

% Now let's visualize the map by plotting the image of
% a square grid of 12 vertical and 6 horizontal lines.

pause     % Strike any key to begin plot

hpplot(w,beta,x,c,12,6)

% Note how the lines intersect at right angles.  Also, 
% note how they converge at the last vertex, the origin,
% since that is the image of the point at infinity.

pause     % Strike any key to continue
clc
% Let's change the fundamental domain from the half-plane 
% to the unit disk.

[z,c]=hp2disk(w,beta,x,c);

ddisp(w,beta,z,c)

pause     % Strike any key to continue

% What's the inverse image of the point -.3-.3i?

zp = dinvmap(-.3-.3i,w,beta,z,c)

pause     % Strike any key to continue

% We should get at least 8 accurate digits, by default.

abs(-.3-.3i - dmap(zp,w,beta,z,c))

pause     % Strike any key to continue

% We can change the map so that -.3-.3i is the image of 0...

[z,c] = dfixwc(w,beta,z,c,-.3-.3i);

pause     % Strike any key to continue

% ...and look at the resulting image of a certain polar grid.

dplot(w,beta,z,c,.1:.1:.9,pi*(.25:.25:2))

pause     % Strike any key to continue
clc
% Now suppose we want an exterior map.  First, the vertices
% have to be given in clockwise order:

w = flipud(w);
beta = scangle(w);

pause     % Strike any key to solve the parameter problem

[z,c] = deparam(w,beta);

pause     % Strike any key to see results

dedisp(w,beta,z,c)

deplot(w,beta,z,c)

echo off      % End of demo

