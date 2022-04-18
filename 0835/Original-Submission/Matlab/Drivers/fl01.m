function [p,z] = fl01()
%
%  generalization of Farmer-Loizou example:
%
%    (x-1)^4k * (x-2)^3k * (x-3)^2K * (X-4)^k
%  
   k = 1;
   p = poly([ones(1,4*k),2*ones(1,3*k),3*ones(1,2*k),4*ones(1,k)]);
   z = [4,k;3,2*k;2,3*k;1,4*k];
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
