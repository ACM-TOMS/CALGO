function [f,variables,zeros,options]=DZ1()
% Generate the DZ1 system, i.e.,
%     x^4_1 - x_2x_3x_4 = 0, 
%     x^4_2 - x_1x_3x_4 = 0, 
%     x^4_3 - x_1x_2x_4 = 0, 
%     x^4_4 - x_1x_2x_3 = 0 at (0, 0, 0, 0).
% The calling sequence is:
%     [f,variables,zeros,options]=DZ1()
% where f, variables, zeros and options can be used by multiplicity directly 
% to compute the multiplicity structure for the DZ1 system. See 
% Test_DZ1 for details.

f{1}='x1^4-x2*x3*x4';
f{2}='x2^4- x1*x3*x4';
f{3}='x3^4-x1*x2*x4';
f{4}='x4^4-x1*x2*x3';
variables{1}='x1';
variables{2}='x2';
variables{3}='x3';
variables{4}='x4';

zeros=[0 0 0 0];
options=optset('EqnType','Poly');
