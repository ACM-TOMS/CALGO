function [p,z] = toh05
%
%  generalization of K.C. Toh and L. N. Trefethen
% 
   k = [-19:0]';
   z = 2.^k;
   p = poly(z);
   z = [z,ones(20,1)];

   if norm(imag(z(:,1))) == 0 
        fprintf('    Illconditioned polynomial, constructed with\n');
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
