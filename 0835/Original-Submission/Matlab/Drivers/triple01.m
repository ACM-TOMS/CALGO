function [p,z] = triple01
%
%  test polynomial suggested by Goedecker
%
    p = triple(5,5,5);
    z = [[0.9,1,1.1]',[5,5,5]'];
    fprintf('                 roots         multiplicities\n');
        fprintf('\n');

        fprintf('%25.15f \t \t \t %3g \n', z');
