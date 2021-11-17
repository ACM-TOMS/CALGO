function [p,z] = petk01
%
% M. Petkovic testing polynomials (p. 109)
%
    p1 = poly([-1, -1, 3, 3, 3, -i*[1,1,1,1]]);
    p2 = [1,-2,5];
    p2 = conv(p2,p2);
    p = conv(p1,p2);
    z = [-1,2; 3, 3; -i, 4; 1+2*i, 2; 1-2*i, 2];

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
