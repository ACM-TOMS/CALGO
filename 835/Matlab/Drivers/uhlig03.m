function [p,z] = uhlig03
%
%  test polynomial used by F. Uhlig
%
   r = [(3/11)*ones(1,12),11/3,11/3,(2*i/7)*ones(1,4),...
           (2.5+i/4)*ones(1,2),1/4];
   p = poly(r);
   z = [3/11,12; 11/3, 2; 2*i/7, 4; 2.5+i/4, 2; 1/4, 1];
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
