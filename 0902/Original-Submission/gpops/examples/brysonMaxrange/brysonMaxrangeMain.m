% ----------------------------
% Bryson Maximum-Range Example
% ----------------------------
% This problem is taken verbatim from the following reference:
% Bryson, A. E. and Ho, Y-C, Applied Optimal Control, Hemisphere
% Publishing, New York, 1975.
% ----------------------------
clear setup limits guess;

global constants
constants.gravity = 1;
constants.acc = 0.5*constants.gravity;

x0 = 0;
y0 = 0;
v0 = 0;
yf = 0.1;
xmin = -10;
xmax = -xmin;
ymin = xmin;
ymax = xmax;
vmin = -100;
vmax = -vmin;
u1min = -10;
u1max = -u1min;
u2min = u1min;
u2max = u1max;
t0 = 0;
tf = 2;

iphase = 1;
limits(iphase).nodes = 40;
limits(iphase).time.min = [t0 tf];
limits(iphase).time.max = [t0 tf];
limits(iphase).state.min(1,:)   = [x0 xmin xmin];
limits(iphase).state.max(1,:)   = [x0 xmax xmax];
limits(iphase).state.min(2,:)   = [y0 ymin yf];
limits(iphase).state.max(2,:)   = [y0 ymax yf];
limits(iphase).state.min(3,:)   = [v0 vmin vmin];
limits(iphase).state.max(3,:)   = [v0 vmax vmax];
limits(iphase).control.min(1,:) = u1min;
limits(iphase).control.max(1,:) = u1max;
limits(iphase).control.min(2,:) = u2min;
limits(iphase).control.max(2,:) = u2max;
limits(iphase).parameter.min    = [];
limits(iphase).parameter.max    = [];
limits(iphase).path.min    = 1;
limits(iphase).path.max    = 1;
limits(iphase).event.min   = [];
limits(iphase).event.max   = [];
guess(iphase).time =  [t0; tf];
guess(iphase).state(:,1) = [x0; x0];
guess(iphase).state(:,2) = [y0; yf];
guess(iphase).state(:,3) = [v0; v0];
guess(iphase).control(:,1) = [0; 1];
guess(iphase).control(:,2) = [1; 0];
guess(iphase).parameter = [];

setup.name = 'Bryson-Maxrange';
setup.funcs.cost = 'brysonMaxrangeCost';
setup.funcs.dae = 'brysonMaxrangeDae';
setup.limits = limits;
setup.derivatives = 'automatic';
setup.checkDerivatives = 1;
setup.guess = guess;
setup.linkages = [];
setup.autoscale = 'off';
output = gpops(setup);
solution = output.solution;
