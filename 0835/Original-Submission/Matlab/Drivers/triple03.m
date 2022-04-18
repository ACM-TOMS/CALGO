function [p,z] = triple03
%
%  test polynomial suggested by Goedecker
%
    p = triple(18,10,16);
    z = [[0.9,1,1.1]',[18,10,16]'];
    fprintf('                 roots         multiplicities\n');
        fprintf('\n');
        fprintf('%25.15f \t \t \t %3g \n', z');

        fprintf('   multroot works when tol = 1.0d-9');
