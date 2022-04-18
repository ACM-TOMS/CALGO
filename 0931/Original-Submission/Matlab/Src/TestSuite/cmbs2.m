function [f,variables,zeros,options]=cmbs2()
%  Generate the cmbs2 system, i.e.,
%    x^3 - 3x^2y + 3xy^2 - y^3 - z^2=0, 
%    z^3 - 3z^2x + 3zx^2 - x^3 - y^2=0, 
%    y^3 - 3y^2z + 3yz^2 - z^3 - x^2=0  at (0, 0, 0).
% The calling sequence is:
%     [f,variables,zeros,options]=cmbs2()
% where f, variables, zeros and options can be used by multiplicity directly 
% to compute the multiplicity structure for the cmbs2 system. See Test_cmbs2 
% for details.

f{1}='x^3 -3*x^2*y + 3*x*y^2 -y^3 - z^2';
f{2}='z^3 - 3*z^2*x + 3*z*x^2 - x^3 - y^2';
f{3}='y^3 - 3*y^2*z + 3*y*z^2 - z^3 - x^2';

variables{1}='x';
variables{2}='y';
variables{3}='z';
zeros=[0 0 0];
options=optset('EqnType','Poly');
