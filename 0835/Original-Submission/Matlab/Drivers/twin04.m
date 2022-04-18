function [p,z] = twin04()

   k = 16;
   p = poly([-0.2*ones(1,k),0.39*ones(1,k),0.40*ones(1,k)]);
   z = [-0.2, k; 0.39, k; 0.4, k];
   fprintf('                 roots         multiplicities\n');
        fprintf('\n');

        fprintf('%25.15f \t \t \t %3g \n', z');
