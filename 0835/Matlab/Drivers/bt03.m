function [p,z] = bt03
%
% Brugnano and Trigiante
%
    p = poly([i*ones(1,5),-i*ones(1,5),0.5i*ones(1,4),...
            -0.5i*ones(1,4),0.75i,-0.75i]);
    z = [i,5; -i, 5; 0.5i, 4; -0.5i, 4; .75i, 1; -0.75i, 1];
    
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
