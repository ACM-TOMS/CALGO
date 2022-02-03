function [f,variables,zeros,options]=mth191()
% Generate the mth191 system, i.e.,
%      x^3 + y^2 + z^2 - 1 = 0, 
%      x^2 + y^3 + z^2 - 1 = 0, 
%      x^2 + y^2 + z^3 - 1 = 0 at (0, 1, 0).
% The calling sequence is:
%     [f,variables,zeros,options]=mth191()
% where f, variables, zeros and options can be used by multiplicity directly 
% to compute the multiplicity structure for the mth191 system. See 
% Test_mth191 for details.

f{1}='x^3 + y^2+z^2-1';
f{2}='y^3 +x^2+z^2-1';
f{3}='z^3 + x^2+y^2-1';
variables{1}='x';
variables{2}='y';
variables{3}='z';
zeros=[0 1 0];
options=optset('EqnType','Poly');
