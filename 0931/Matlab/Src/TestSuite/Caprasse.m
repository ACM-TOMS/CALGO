function [f,variables,zeros,options]=Caprasse()
% Generate the Caprasse system.  The calling sequence is:
%     [f,variables,zeros,options]=Caprasse()
% where f, variables, zeros and options can be used by multiplicity directly 
% to compute the multiplicity structure for the Caprasse system. See 
% Test_Caprasse for details.

f{1}='-x1^3*x3+4*x1*x2^2*x3+4*x1^2*x2*x4+2*x2^3*x4+4*x1^2-10*x2^2+4*x1*x3-10*x2*x4+2';
f{2}='-x1*x3^3 + 4*x2*x3^2*x4 + 4*x1*x3*x4^2 + 2*x2*x4^3 + 4*x1*x3 + 4*x3^2- 10*x2*x4-10*x4^2 + 2';
f{3}='x2^2*x3+2*x1*x2*x4-2*x1-x3';
f{4}='2*x2*x3*x4+x1*x4^2-x1-2*x3';
variables{1}='x1';
variables{2}='x2';
variables{3}='x3';
variables{4}='x4';

zeros=[2 -1i*sqrt(3) 2 1i*sqrt(3)];
options=optset('EqnType','Poly');
