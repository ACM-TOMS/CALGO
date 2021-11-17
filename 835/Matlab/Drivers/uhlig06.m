function [p,z] = uhlig06
%
%  test polynomial used by F. Uhlig
%
   p = [1,0,0,0,0,0,0,0,-1];
   z = [-1.00000000000000                    
 -0.70710678118655 + 0.70710678118655i
 -0.70710678118655 - 0.70710678118655i
                 0 + 1.00000000000000i
                 0 - 1.00000000000000i
  1.00000000000000                    
  0.70710678118655 + 0.70710678118655i
  0.70710678118655 - 0.70710678118655i];
    z = [z,ones(8,1)];
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
