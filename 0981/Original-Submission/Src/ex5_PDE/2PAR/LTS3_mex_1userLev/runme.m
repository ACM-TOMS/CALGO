% runme.m
%{
	*** MATLAB script file to build and run a sample program using
	***              Talbot Suite DE
	***
	*** author: Mariarosaria Rizzardi
	***         mariarosaria.rizzardi@uniparthenope.it
	***
%}
clear; clc
choice = menu('Select:','1)  accuracy test','2)  time test');

switch choice
    case 1 % ACCURACY TEST
        fprintf('entered 1: accuracy test\n')
        fprintf('\nBUILD EXECUTABLE FOR ACCURACY TEST\n\tbuild_accuracy\n\n')
        build_accuracy
        fprintf('\n********************************************\n\n')
        tol=1e-12; thrds1=4; thrds2=2;
        fprintf('RUN EXECUTABLE FOR ACCURACY TEST\n\ttol=%6.0e; thrds1=%d; thrds2=%d; OMP_main_ACCURACY(tol,thrds1,thrds2)\n', tol,thrds1,thrds2)
        OMP_main_ACCURACY(tol,thrds1,thrds2)
        fprintf('\n\n********************************************\n')
        fprintf('\tTO REMOVE EXECUTABLE FILES,\n')
        fprintf('\tENTER:\n')
        fprintf('\t\tclean\n')
        fprintf('\tOTHERWISE RUN AS:\n\t\ttol=%6.0e; thrds1=%d; thrds2=%d; OMP_main_ACCURACY(tol,thrds1,thrds2)', tol,thrds1,thrds2)
        fprintf('\n\n********************************************\n\n')
    otherwise % TIME TEST
        fprintf('entered 2: time test\n')
        fprintf('\nBUILD EXECUTABLE FOR TIME TEST\n\tbuild_times\n\n')
        build_times
        fprintf('\n********************************************\n\n')
        tol=1e-12; jFUN=1; NTval=120; NXYval=120;
        fprintf('RUN EXECUTABLE FOR TIME TEST\n\ttol=%6.0e; jFUN=%d; NTval=%d; NXYval=%d; OMP_main_TIMES(tol,jFUN,NTval,NXYval)\n', tol,jFUN,NTval,NXYval)
        OMP_main_TIMES(tol,jFUN,NTval,NXYval)
        fprintf('\n\n********************************************\n')
        fprintf('\tTO REMOVE EXECUTABLE FILES,\n')
        fprintf('\tENTER:\n')
        fprintf('\t\tclean\n')
        fprintf('\tOTHERWISE RUN AS:\n\t\ttol=%6.0e; jFUN=%d; NTval=%d; NXYval=%d; OMP_main_TIMES(tol,jFUN,NTval,NXYval)', tol,jFUN,NTval,NXYval)
        fprintf('\n\n********************************************\n\n')
end
