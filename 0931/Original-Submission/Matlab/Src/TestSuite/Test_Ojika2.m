function [m,DI,HF]=Test_Ojika2(flag)
% Compute the multiplicity structure for Ojika2 system, i.e.,
%     x^2 + y + z - 1 = 0, 
%     x + y^2 + z - 1 = 0, 
%     x + y + z^2 - 1 = 0 at (0, 0, 1) or (1, 0, 0).
% The calling sequence is:
%     [m,DI,HF]=Test_Ojika2(flag)
%   flag indicates the zero to be tested, i.e,
%        1 for (0,0,1) and 2 for (1,0,0).
%Output:
%   m returns the multiplicity.
%   DI returns the basis for the dual space.
%   HF returns the Hilbert function.

[f,variables,zeros,options]=Ojika2(flag);
disp('====================Ojika2 system===================');
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
disp(' ')
pause(1)
disp(['   multiplicity: ' num2str(m)])
pause(1)
disp('   the dual basis: ')
outputDI(DI)
pause(1)
disp(['   the Hilbert function: ' num2str(HF)])
disp(' ')

