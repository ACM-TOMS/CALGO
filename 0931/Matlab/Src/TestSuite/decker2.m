function [f,variables,zeros,options]=decker2()
% Generate the decker2 system, i.e.,
%      x + y^3=0, x^2y - y^4=0, at (0, 0).
% The calling sequence is:
%     [f,variables,zeros,options]=decker2()
% where f, variables, zeros and options can be used by multiplicity directly 
% to compute the multiplicity structure for the decker2 system. See 
% Test_decker2 for details.

f=[sym('x + y^3'); sym('x^2*y - y^4')];
variables=[sym('x');sym('y')];
zeros=[0 0];
options=optset('EqnType','Poly');
