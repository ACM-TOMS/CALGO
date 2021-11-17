function [m,DI,HF]=Test_decker2()
% Compute the multiplicity structure for decker2 system, i.e.,
%      x + y^3=0, x^2y - y^4=0, at (0, 0).
% The calling sequence is:
%     [m,DI,HF]=Test_decker2()
%
%Output:
%   m returns the multiplicity.
%   DI returns the basis for the dual space.
%   HF returns the Hilbert function.

[f,variables,zeros,options]=decker2();
disp('====================decker2 system===================');
disp(['   system: ' char(f(1))])
for n=2:length(f)
    disp(['           ' char(f(n))])
end
pause(1)
disp(['   zero: ' num2str(zeros)])
pause(1)
disp(' ')
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

