function [f,variables,zeros,options]=cmbs1()
% Generate the cmbs1 system, i.e.,
%    x^3 - yz = 0,
%    y^3 - xz = 0, 
%    z^3 - xy = 0 at (0, 0, 0).
% The calling sequence is:
%     [f,variables,zeros,options]=cmbs1()
% where f, variables, zeros and options can be used by multiplicity directly 
% to compute the multiplicity structure for the cmbs1 system. See Test_cmbs1 
% for details.

f{1}='x^3 - y*z';
f{2}='y^3 - x*z';
f{3}='z^3 - x*y';
variables{1}='x';
variables{2}='y';
variables{3}='z';
zeros=[0 0 0];
options=optset('EqnType','Poly');