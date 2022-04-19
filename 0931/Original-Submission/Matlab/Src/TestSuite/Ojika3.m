function [f,variables,zeros,options]=Ojika3(flag)
% Generate the Ojika3 system, i.e.,
%      x + y + z - 1 = 0, 
%      2x^3 + 5y^2 - 10z + 5z^3 + 5 = 0, 
%      2x + 2y + z^2 - 1 = 0 at (0, 0, 1) or (-5/2,5/2, 1).
% The calling sequence is:
%     [f,variables,zeros,options]=Ojika3(flag)
% where flag is the type of zero (1 for (0, 0, 1) and 2 for (-5/2,5/2, 1)), 
% f, variables, zeros and options can be used by multiplicity directly to 
% compute the multiplicity structure for the Ojika3 system. See Test_Ojika3 
% for details.

f{1}='x +y+z-1';
f{2}='2*x^3 + 5*y^2 - 10*z + 5*z^3 + 5';
f{3}='2*x + 2*y + z^2 - 1';
variables{1}='x';
variables{2}='y';
variables{3}='z';
if flag==1
    zeros=[0 0 1];
elseif flag==2
    zeros=[-5/2 5/2 1];
end
options=optset('EqnType','Poly');
