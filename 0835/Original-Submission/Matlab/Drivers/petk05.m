function [p,z] = petk05
%
% M. Petkovic testing polynomials, p142
%
    p = [1,-7,20,-28,-18,110,-92,-44,345,225];
    z = [3.00000000000000 + 0.00000000000000i
  1.00000000000000 + 2.00000000000000i
  1.00000000000000 - 2.00000000000000i
 -1.00000000000000 - 0.00000000000000i];
    z = [z, [2,2,2,3]'];
    

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
