%% Two dielectric cylinders (obstacle to Fresnel in 2D: twodielTM_8f.exp)
% For reference see <incontrastsRef.html>.
%
%% Positions of the obstacles
%
% * Theoretical, see [1], position should be left and right from origin, but Fresnel 
%   data are rotated 270 degrees mathematical positive.
% * Actually the cylinders are above and below the origin.
%
%% Other parameters of the obstacles
% See <fresnel_op1_dielTM.html> (because it is a dielectric too).

%% References
%
% * [1] Kamal Belkebir and Marc Saillard. Special section on testing inversion algorithms against experimental data. _Inverse Problems_, 17(6):1565-1571, 2001.
% * [2] Florian B&uuml;rgel, Kamil S. Kazimierski, and Armin Lechleiter. A sparsity regularization and total variation based computational framework for the inverse medium problem in scattering. _Journal of Computational Physics_, 339:1-30, 2017.

function q = fresnel_op1_twodielTM(X1,X2,varargin)
% two dielectric cylinders

qBall = 2 + 1i*0; % dielectricum, so imaginary part is 0
rBall = 15E-3; % 15 mm

% correct the position:

positionmethod = 1; % 1 or 2: both are manually corrections of the obstacles position
% method 1 is used in [2].

switch positionmethod

  case 1 % method 1--------------------------------------------------
% 1) a little shift by shifting the grid

xshift = -8E-3;
yshift = +1E-3;

X1 = X1-xshift;
X2 = X2-yshift;

% 2) rotation by rotation of grid in this file:

gamma = 8; %rotation angle in degrees (mathematical positive)

%rotation matrix:
R = @(alpha) [cos(alpha), -sin(alpha); sin(alpha), cos(alpha)]; % alpha is angle in rad...

Xrot = R(gamma*2*pi/360)*[X1;X2];
X1 = Xrot(1,:);
X2 = Xrot(2,:);

%-----------------------------------------------------

dist = 45E-3; % distance from origin

r1 = hypot(X1,X2-dist); % 45 mm below origin
r2 = hypot(X1,X2+dist); % 45 mm above origin

% this is changed because of shift and rotation of X1 and X2

  case 2 % method 2------------------------------------------

% set positions manually

p = @(n) -0.06+n/90*0.12;

pos1 = [p(44); p(78)];
pos2 = [p(36); p(15)];

r1 = hypot(X1-pos1(1),X2-pos1(2)); % pos1
r2 = hypot(X1-pos2(1),X2-pos2(2)); % pos2

end

q1 = ( r1 < rBall );
q2 = ( r2 < rBall );
q = qBall*(q1 | q2);
q = double(q);

% Be careful: this does not work... you end up with 0 and 1 entries:
% q1 = qBall*( r1 < rBall ); % OK
% q2 = qBall*( r2 < rBall ); % OK
% q = q1 | q2;               % this resuts in 0 or 1 entries...

end
