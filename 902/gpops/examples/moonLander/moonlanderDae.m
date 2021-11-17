function daeout = moonlanderDae(soldae);

global CONSTANTS

t = soldae.time;
x = soldae.state;
u = soldae.control;
p = soldae.parameter;

h = x(:,1);
v = x(:,2);

hdot = v;
vdot = -CONSTANTS.g+u;

daeout = [hdot vdot];
