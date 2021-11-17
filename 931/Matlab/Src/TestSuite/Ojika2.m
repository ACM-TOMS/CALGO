function [f,variables,zeros,options]=Ojika2(flag)
% Generate the Ojika2 system, i.e.,
%     x^2 + y + z - 1 = 0, 
%     x + y^2 + z - 1 = 0, 
%     x + y + z^2 - 1 = 0 at (0, 0, 1) or (1, 0, 0).
% The calling sequence is:
%     [f,variables,zeros,options]=Ojika2(flag)
% where flag is the type of zero (1 for (0, 0, 1) and 2 for (1, 0, 0)), f, 
% variables, zeros and options can be used by multiplicity directly to 
% compute the multiplicity structure for the Ojika2 system. See Test_Ojika2 
% for details.

f{1}='x^2 +y+z-1';
f{2}='y^2 + x+z-1';
f{3}='z^2 + x+y-1';
variables{1}='x';
variables{2}='y';
variables{3}='z';
if flag ==1
    zeros=[0 0 1];
elseif flag==2
    zeros=[1 0 0];
end
options=optset('EqnType','Poly');
