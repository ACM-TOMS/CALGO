function [m,DI,HF]=Test_Ojika3(flag)
% Compute the multiplicity structure for Ojika3 system, i.e.,
%      x + y + z - 1 = 0, 
%      2x^3 + 5y^2 - 10z + 5z^3 + 5 = 0, 
%      2x + 2y + z^2 - 1 = 0 at (0, 0, 1) or (-5/2,5/2, 1).
% The calling sequence is:
%     [m,DI,HF]=Test_Ojika3(flag)
%   flag indicates the zero to be tested, i.e,
%        1 for (0,0,1) and 2 for (-5/2,5/2,1).
%Output:
%   m returns the multiplicity.
%   DI returns the basis for the dual space.
%   HF returns the Hilbert function.

[f,variables,zeros,options]=Ojika3(flag);
disp('====================Ojika3 system===================');
disp(['   system: ' f{1}])
for n=2:length(f)
    disp(['           ' f{n}])
end
pause(1)
disp(['   zero: ' num2str(zeros)])
pause(1)
disp(' ');
disp('   Start computing the multiplicity structure')
t1=cputime;
[m,DI,HF]=multiplicity(f,variables,zeros,options);
t2=cputime;
disp(['   computation finished in ',num2str(t2-t1),' seconds, with following results:']);
disp(' ');
pause(1)
disp(['   multiplicity: ' num2str(m)])
pause(1)
disp('   the dual basis: ')
outputDI(DI)
pause(1)
disp(['   the Hilbert function: ' num2str(HF)])
disp(' ')

