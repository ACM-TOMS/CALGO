function [p,z] = near03
%
%  test polynomial suggested by Z. Zeng
%  It is suggested that MultRoot input gamma = 10.
%
    e = 0.001;
    p = poly([(1-e)*ones(1,20),ones(1,20),-0.5*[1,1,1,1,1]]);
    z = [1-e, 20; 1, 20; -0.5, 5];
    if norm(imag(z(:,1))) == 0 
        fprintf('                 roots         multiplicities\n');
        fprintf('\n');
        fprintf('%25.15f \t \t \t %3g \n', z');
    else
        fprintf('                 roots ')
        fprintf('   \t\t\t\t\t\t     multiplicities\n');
        fprintf('\n');
        fprintf('%22.15f + %22.15f i \t \t %3g \n', ...
            [real(z(:,1)),imag(z(:,1)),z(:,2)]');
    end;    
    fprintf('To use MultRoot, set input gamma = 5\n')
