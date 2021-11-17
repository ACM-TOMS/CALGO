function TestAll

disp('SMOOTH KERNEL')
TestSmooth(8,1,-pi,pi)
pause
disp(' ')
disp('********************************************************************')
disp(' ')

disp('DISCONTINUOUS KERNEL')
TestDiscontinuous(3,1,0,5)
pause
disp(' ')
disp('********************************************************************')
disp(' ')

disp('LOGARITHMIC SINGULARITY')
TestLogSing(5,1,0,2)
pause
disp(' ')
disp('********************************************************************')
disp(' ')

disp('ALGEBRAIC SINGULARITY')
TestAlgSing(2,5,-1,1,0.5)
pause
disp(' ')
disp('********************************************************************')
disp(' ')
disp('INFINITE INTERVAL')
TestInf(14,2,[],[],15)