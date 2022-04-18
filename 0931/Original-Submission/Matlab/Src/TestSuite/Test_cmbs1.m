function [m,DI,HF]=Test_cmbs1()
% Compute the multiplicity structure for cmbs1 system, i.e.,
%    x^3 - yz = 0,
%    y^3 - xz = 0, 
%    z^3 - xy = 0 at (0, 0, 0).
% The calling sequence is:
%     [m,DI,HF]=Test_cmbs1()
%
%Output:
%   m returns the multiplicity.
%   DI returns the basis for the dual space.
%   HF returns the Hilbert function.

[f,variables,zeros,options]=cmbs1();
disp('====================cmbs1 system===================');
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
