function [p,z] = toh06c
%
%  generalization of K.C. Toh and L. N. Trefethen
% 
   p = ones(1,6);
   p = conv(p,p); p = conv(p,p);
     z = [    0.50000000000000 + 0.86602540378444i
  0.50000000000000 - 0.86602540378444i
 -1.00000000000000                    
 -0.50000000000000 + 0.86602540378444i
 -0.50000000000000 - 0.86602540378444i];
    z = [z,4*ones(5,1)];
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
