function [p,z] = triple02
%
%  test polynomial suggested by Goedecker
%
    p = triple(10,10,10);
    z = [[0.9,1,1.1]',[10,10,10]'];
    fprintf('                 roots         multiplicities\n');
        fprintf('\n');

        fprintf('%25.15f \t \t \t %3g \n', z');
