% Schwarz-Christoffel Toolbox
% Version 1.3   June 1, 1995.
% Written by Toby Driscoll (driscoll@na-net.ornl.gov).
% See user's guide for full usage details.
%
% Graphical user interface (GUI).
%   scgui      - Activate graphical user interface.
%   scgget     - Get polygon/solution properties from GUI.
%   scgset     - Set polygon/solution properties in the GUI.
% 
% Working with polygons.
%   drawpoly   - Draw a polygon with the mouse.
%   plotpoly   - Plot a polygon.
%   modpoly    - Modify a polygon with the mouse.
%   scselect   - Select polygon vertices with the mouse.
%   scaddvtx   - Add a vertex.
%   scangle    - Compute turning angles.
%   scfix      - Make polygon acceptable to other SC routines.
% 
% Half-plane->polygon map.
%   hpparam    - Solve parameter problem.
%   hpdisp     - Pretty-print solution data.
%   hpmap      - Compute forward map.
%   hpinvmap   - Compute inverse map.
%   hpplot     - Adaptive plotting of the image of a cartesian grid.
%   hpderiv    - Derivative of the map.
% 
% Disk->polygon map.
%   dparam     - Solve parameter problem.
%   ddisp      - Pretty-print solution data.
%   dmap       - Compute forward map.
%   dinvmap    - Compute inverse map.
%   dplot      - Adaptive plotting of the image of a polar grid.
%   dderiv     - Derivative of the map.
% 
% Disk->exterior polygon map.
%   deparam    - Solve parameter problem.
%   dedisp     - Pretty-print solution data.
%   demap      - Compute forward map.
%   deinvmap   - Compute inverse map.
%   deplot     - Adaptive plotting of the image of a polar grid.
%   dederiv    - Derivative of the map.
% 
% Strip->polygon map.
%   stparam    - Solve parameter problem.
%   stdisp     - Pretty-print solution data.
%   stmap      - Compute forward map.
%   stinvmap   - Compute inverse map.
%   stplot     - Adaptive plotting of the image of a polar grid.
%   stderiv    - Derivative of the map.
% 
% Rectangle->polygon map.
%   rparam    - Solve parameter problem.
%   rdisp     - Pretty-print solution data.
%   rmap      - Compute forward map.
%   rinvmap   - Compute inverse map.
%   rplot     - Adaptive plotting of the image of a polar grid.
%   rderiv    - Derivative of the map.
% 
% Conversion routines.
%   hp2disk    - Convert a solution from half-plane to one from disk.
%   disk2hp    - Convert a solution from disk to one from half-plane.
%   dfixwc     - Choose conformal center of disk map.
%   ptsource   - Graphical use of DFIXWC.
%
% Demonstrations.
%   scdemo     - Select demos from a menu.
%   tutdemo    - Walk through a tutorial.
%   infdemo    - Explain infinite vertices.
%   elongdemo  - Maps to elongated polygons.
%   faberdemo  - Introduce Faber polynomials.
