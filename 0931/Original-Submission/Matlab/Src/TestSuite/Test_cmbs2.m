function [m,DI,HF]=Test_cmbs2()
% Compute the multiplicity structure for cmbs2 system, i.e.,
%    x^3 - 3x^2y + 3xy^2 - y^3 - z^2=0, 
%    z^3 - 3z^2x + 3zx^2 - x^3 - y^2=0, 
%    y^3 - 3y^2z + 3yz^2 - z^3 - x^2=0  at (0, 0, 0).
% The calling sequence is:
%     [m,DI,HF]=Test_cmbs2()
%
%Output:
%   m returns the multiplicity.
%   DI returns the basis for the dual space.
%   HF returns the Hilbert function.

[f,variables,zeros,options]=cmbs2();
disp('====================cmbs2 system===================');
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

