function [p,z] = petk02
%
% M. Petkovic testing polynomials, p118
%
    p = poly([1,1,-i,-i,-i,5*i,5*i,-5*i,-5*i]);
    z = [1,2; -i,3; 5*i,2; -5*i,2];

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
