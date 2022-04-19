function [p,z] = jt10b
%
%  test polynomial suggested by Jenkins and Traub
%
    a = 10^6;
    p = poly([a,1,1/a]);
    z = [[a,1,1/a]',ones(3,1)];
    fprintf('\n');
    if norm(imag(z(:,1))) == 0 
        fprintf('                 roots         multiplicities\n');
        fprintf('\n');
        fprintf('%25.8f \t \t \t %3g \n', z');
    else
        fprintf('                 roots ')
        fprintf('   \t\t\t\t\t\t     multiplicities\n');
        fprintf('\n');
        fprintf('%22.15f + %22.15f i \t \t %3g \n', ...
            [real(z(:,1)),imag(z(:,1)),z(:,2)]');
    end;     
