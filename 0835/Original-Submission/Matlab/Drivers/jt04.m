function [p,z] = jt04
%
%  test polynomial suggested by Jenkins and Traub
%
   p = poly([0.1,0.1,0.1,0.5,0.6,0.7]);
   z = [0.5, 1; 0.6, 1; 0.7, 1; 0.1, 3];
   fprintf('\n');
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
