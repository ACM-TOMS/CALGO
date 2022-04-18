function [p,z] = toh06a
%
%  generalization of K.C. Toh and L. N. Trefethen
% 
   p = ones(1,11);
   p = conv(p,p);
     z = [  0.84125353283118 + 0.54064081745560i
  0.84125353283118 - 0.54064081745560i
  0.41541501300189 + 0.90963199535452i
  0.41541501300189 - 0.90963199535452i
 -0.95949297361449 + 0.28173255684143i
 -0.95949297361449 - 0.28173255684143i
 -0.65486073394529 + 0.75574957435426i
 -0.65486073394529 - 0.75574957435426i
 -0.14231483827329 + 0.98982144188093i
 -0.14231483827329 - 0.98982144188093i];
    z = [z,2*ones(10,1)];
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
