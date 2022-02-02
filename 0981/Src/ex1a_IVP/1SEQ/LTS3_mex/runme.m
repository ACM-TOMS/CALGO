% runme.m
%{
	*** MATLAB script file to build and run a sample program using
	***              Talbot Suite DE
	***              Example 1a - LT samples by ode45.m
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
        fprintf('RUN EXECUTABLE FOR ACCURACY TEST\n\ttol=1e-6; SEQ_main_ACCURACY(tol)\n')
        tol=1e-12; SEQ_main_ACCURACY(tol)
        fprintf('\n\n********************************************\n')
        fprintf('\tTO REMOVE EXECUTABLE FILES,\n')
        fprintf('\tENTER:\n')
        fprintf('\t\tclean\n')
        fprintf('\tOTHERWISE RUN AS:\n')
        fprintf('\t\ttol=1e-6; SEQ_main_ACCURACY(tol)')
        fprintf('\n\n********************************************\n\n')
    otherwise % TIME TEST
        fprintf('entered 2: time test\n')
        fprintf('\nBUILD EXECUTABLE FOR TIME TEST\n\tbuild_times\n\n')
        build_times
        fprintf('\n********************************************\n\n')
        fprintf('RUN EXECUTABLE FOR TIME TEST\n\ttol=1e-6; SEQ_main_TIMES(tol)\n')
        tol=1e-12; SEQ_main_TIMES(tol)
        fprintf('\n\n********************************************\n')
        fprintf('\tTO REMOVE EXECUTABLE FILES,\n')
        fprintf('\tENTER:\n')
        fprintf('\t\tclean\n')
        fprintf('\tOTHERWISE RUN AS:\n')
        fprintf('\t\ttol=1e-6; SEQ_main_TIMES(tol)')
        fprintf('\n\n********************************************\n\n')
end
