function [p,z] = jt05
%
%  test polynomial suggested by Jenkins and Traub
%
   p = poly([0.1*ones(1,4),0.2*ones(1,3),0.3,0.3,0.4]);
   z = [0.4, 1; 0.3, 2; 0.2, 3; 0.1, 4];
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
