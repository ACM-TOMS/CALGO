function [p,z] = triple04
%
%  test polynomial suggested by Goedecker
%
    p = triple(20,15,10);
    z = [[0.9,1,1.1]',[20,15,10]'];
    fprintf('                 roots         multiplicities\n');
        fprintf('\n');
        fprintf('%25.15f \t \t \t %3g \n', z');

        fprintf('   multroot should work when tol = 1.0d-12');
