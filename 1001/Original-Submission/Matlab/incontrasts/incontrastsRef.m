%% Predefined contrasts in folder "incontrasts"
% Reference to the predefined contrasts in folder "incontrasts" and their 
% usage in the computational framework.
%
%
%% Syntax
% Example of syntax for a ball in 2D and 3D. You can find a list of all 
% predefined contrasts below.
%
%   q = referenceBall2D(X1,X2,seti)     % 2D
%   q = referenceBall3D(X1,X2,X3,seti)  % 3D
%
%
%% Description
% |q = referenceBall2D(X1,X2,seti)| evaluates the contrast in the full grid
% described by X1 and X2 and stores it as a vector 1 x n. Note that X1 and X2 
% are produced with meshgrid, but stored as vectors 1 x n.
%
% |q = referenceBall3D(X1,X2,X3,seti)| evaluates the contrast in the full grid
% described by X1, X2 and X3 and stores it as a vector 1 x n. Note that X1, X2, and X3 
% are produced with meshgrid, but stored as vectors 1 x n.
%
%
%% Input Arguments
%
% * X1, X2, (X3)    : coordinate arrays X1, X2, and X3 (last one in 3D case)
%                     (vectors of size 1 x seti.nROI^seti.dim)
% * seti            : struct for additional input, e.g.
%                     contrast value |seti.qBall| for ball,
%                     radius |seti.rBall| for ball.
%
%
%% Output Arguments
% * q               : contrast evaluated on the region of interet (ROI) 
%                     stored as vector 1 x n.
%
%% Examples in 2D
%
% *Triangle in 2D*
%
% For this example make sure to be in the folder code/incontrasts/2D.
%
%   % Prepare a grid
%   h = 0.016;
%   r = 0.2;
%   x1 = (-r:h:r-h);
%   [X1,X2] = meshgrid(x1,x1);
%   X1 = transpose(X1(:));
%   X2 = transpose(X2(:));
%   
%   % Evaluate predefined contrast corner2D on the grid an plot result
%   seti.rCD = 2*r;
%   q = triangle2D(X1,X2,seti);
%   l = length(x1);
%   qMatrix = reshape(q,[l l]);
%   imagesc(real(qMatrix)); axis xy; colorbar;
%
%
% *Ball in 2D*
%
% For this example make sure to be in the folder code.
%
%   % add all pathes of rebis (needed in case of checkfield)
%   init;
%
%   % Prepare a grid
%   h = 0.016;
%   r = 0.2;
%   x1 = (-r:h:r-h); 
%   [X1,X2] = meshgrid(x1,x1);
%   X1 = transpose(X1(:));
%   X2 = transpose(X2(:));
%
%   % set strcut seti
%   seti.qBall = 0.8;
%   seti.rBall = 0.1;
%
%   % Evaluate predefined contrast referenceBall2D on the grid an plot result
%   q = referenceBall2D(X1,X2,seti);
%   l = length(x1);
%   qMatrix = reshape(q,[l l]);
%   imagesc(real(qMatrix)); axis xy; colorbar;
%
%
%% List of all predefined Contrasts
%
% * q is the contrast value of the described obstacle.
% * Imaginary part of contrast is in all predefined contrasts is 0, if nothing
% other is mentioned.
%
% Folders in incontrasts:
%
% * 2D: several two dimensional predefined contrasts.
% * 2DFresnel: two dimensional predefined contrasts fitting to some of the
% real-world data from Institute Fresnel.
% * 3D: several three dimensional predefined contrasts.
%
% *Predefined contrasts in 2D*
%
% Folder: incontrasts/2D
%
% <html>
% <table>
% <tr><td><em>String</em></td><td><em>Description</em></td></tr>
% <tr><td><a href="empty2D.html">empty2D</a></td><td><strong>Empty contrast</strong> with q = 0 everywhere (to test)</td></tr>
% <tr><td><a href="corner2D.html">corner2D</a></td><td><strong>Corner</strong><br> top and left line with q = 1</td></tr>
% <tr><td><a href="cornerBallSparse2D.html">cornerBallSparse2D</a></td><td><strong>Non-constant corner, ball (filled circle) and broken corner</strong><br>Corner with non-constant contrast q<br> Ball with q = 1 <br>Broken corner with q = 1 (corner consisting of squares and a thin line - very sparse)</td></tr>
% <tr><td><a href="cornerBallSparseMod2D.html">cornerBallSparseMod2D</a></td><td><strong>Non-constant corner, ball (filled circle) and broken corner</strong><br>Same as cornerBallSparse2D, but rotated 180 degrees.</td></tr>
% <tr><td><a href="cross2D.html">cross2D</a></td><td><strong>Cross</strong> with q = 1</td></tr>
% <tr><td><a href="rectangle2D.html">rectangle2D</a></td><td><strong>Rectangle</strong> with q = 0.5</td></tr>
% <tr><td><a href="referenceBall2D.html">referenceBall2D</a></td><td><strong>Ball (filled circle)</strong><br> with contrast value q = seti.qBall (default: 0.8) and radius seti.rBall (default: 0.015)</td></tr>
% <tr><td><a href="referenceBallSmooth2D.html">referenceBallSmooth2D</a></td><td><strong>Smooth ball (filled circle)</strong><br> with contrast value seti.qBall and radius seti.rBall</td></tr>
% <tr><td><a href="shepp2D.html">shepp2D</a></td><td><strong>Modified Shepp-Logan phantom</strong> from MATLAB</td></tr>
% <tr><td><a href="triangle2D.html">triangle2D</a></td><td><strong>Trianlge</strong> with q = 0.8</td></tr>
% <tr><td><a href="twoCorners2D.html">twoCorners2D</a></td><td><strong>Two corners</strong> with q = 1</td></tr>
% <tr><td><a href="twoCornersOneBall2D.html">twoCornersOneBall2D</a></td><td><strong>Two corners and one ball</strong><br>Two corners with q = 0.8<br>Ball with q = 1</td></tr>
% </table>
% </html>
%
% *Predefined contrasts in 2D corresponding to real-world data from Institute Fresnel*
%
% Folder: incontrasts/2DFresnel
%
% The obstacles corresponds to the described ones in [1].
%
% Note that dielectricum means that the imaginary part of the contrast
% vanishes.
%
% <html>
% <table>
% <tr><td><em>String</em></td><td><em>Description</em></td></tr>
% <tr><td><a
% href="fresnel_op1_dielTM.html">fresnel_op1_dielTM</a></td><td>
%   <strong>One dielectric cylinder</strong><br>
%   Contrast value q = 2 + 0 i.<br>
%   Ball with radius 15 mm<br>
%   Distance of center of the ball from origin is approximately 30 mm.<br>
%   (Note:  Position was manually corrected.)
% </td></tr>
% <tr><td><a
% href="fresnel_op1_twodielTM.html">fresnel_op1_twodielTM</a></td><td>
%   <strong>Two dielectric cylinders</strong><br>
%   Contrast value q = 2 + 0 i.<br>
%   Balls with radius 15 mm<br>
%   Distance of center of the balls from origin is approximately 45 mm.<br>
%   (Note:  Position was manually corrected.)
% </td></tr>
% <tr><td><a href="fresnel_op1_rectTM_cent.html">fresnel_op1_rectTM_cent</a></td><td>
%   <strong>Rectangle (experimentally!)</strong><br>
%   This obstacle is experimentally and may contain errors, because the
%   contrast value is assumed because reconstruction.<br>
%   The obstacle is out of metal, i.e. the refractive index is highly
%   dependent on wavelength.<br>
%   Metal is perfectly conducting, i.e. the imaginary part is
%   infinity.<br>
%   In a first approach we assume the contrast value q = 0.5 + 1E6 i.<br>
%   (Note: This contrast is experimentally.)
% </td></tr>
% </table>
% </html>
%
% *Predefined contrasts in 3D*
%
% <html>
% <table>
% <tr><td><em>String</em></td><td><em>Description</em></td></tr>
% <tr><td><a href="empty3D.html">empty3D</a></td><td><strong>Empty contrast</strong> with q = 0 everywhere (to test)</td></tr>
% <tr><td><a href="corner3D.html">corner3D</a></td><td>
%   <strong>Corner (tripod)</strong><br>
%   first horizontal arc with q = 1.0,<br>
%   second horizontal arc with q = 0.8,<br>
%   vertical horizontal arc with q = 0.6.
% </td></tr>
% <tr><td><a href="cross3D.html">cross3D</a></td><td>
%   <strong>Cross</strong><br>
%   first horizontal bar: q = 0.8,<br>
%   second horizontal bar with q = 0.6,<br>
%   vertical bar with q = 1.0.
% </td></tr>
% <tr><td><a href="referenceBall3D.html">referenceBall3D</a></td><td>
%   <strong>Ball</strong><br>
%   with contrast value q = seti.qBall (default: 0.8)
%   and radius seti.rBall (default: 0.015).
% </td></tr>
% <tr><td><a href="twoTripods.html">twoTripods3D</a></td><td>
%   <strong>Two tripods</strong><br>
%   all arcs with q = 1.0.<br>
% </td></tr>
% <tr><td><a href="cubeLike3D.html">cubeLike3D</a></td><td>
%   <strong>Edges of a cube</strong><br>
%   with contrast value q = 1.
% </td></tr>
% </table>
% </html>
%
%
%% Usage in computational framework
% Usage of predefined contrasts in computational framework
%
% *Syntax*
%
%   seti.contrast = str;
%
% *Description*
%
% |seti.contrast = str| uses the contrast defined in function |str|.
%
% * |seti.contrast| is defined in a file in folder |inseti|.
%
% *Example*
%
% Define in inseti/example.m: (see <example.html>)
%
%   seti.contrast = 'cornerBallSparse2D';
%
%
%% References
%
% * [1] Kamal Belkebir and Marc Saillard. Special section on testing inversion algorithms against experimental data. _Inverse Problems_, 17(6):1565-1571, 2001.
%
