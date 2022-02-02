function [p,z] = igyp02c
%
%  generalization of Igarash and Ypma
% 
   m = 3;
   p = poly([10*(1+i)*ones(1,m),-1*ones(1,10-m)]);
   z = [10*(1+i),m; -1, 10-m];
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
