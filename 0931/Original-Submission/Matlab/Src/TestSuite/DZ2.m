function [f,variables,zeros,options]=DZ2()
% Generate the DZ2 system, i.e.,
%           x^4 = 0, 
%           x^2y + y^4 = 0, 
%           z + z^2 - 7x^3 - 8x^2 = 0 at (0, 0, 1)
% The calling sequence is:
%     [f,variables,zeros,options]=DZ2()
% where f, variables, zeros and options can be used by multiplicity directly 
% to compute the multiplicity structure for the DZ2 system. See 
% Test_DZ2 for details.

f{1}='x^4 ';
f{2}='x^2*y + y^4';
f{3}='z+z^2- 7*x^3-8*x^2';
variables{1}='x';
variables{2}='y';
variables{3}='z';
zeros=[0 0 -1];
options=optset('EqnType','Poly');
