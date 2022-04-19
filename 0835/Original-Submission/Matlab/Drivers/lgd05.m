function [p,z] = lgd05
%
%  test polynomial suggested by Goedecker
%   Legendre polynomial of degree 5
%
    p = lgd(5);
    z = [                  0, ...
   0.90617984593866, ...
   0.53846931010568, ...
  -0.90617984593866, ...
  -0.53846931010568];
    z = [z',ones(5,1)];
    fprintf('\n');
   if norm(imag(z(:,1))) == 0 
        fprintf('                 roots         multiplicities\n');
        fprintf('\n');
        fprintf('%25.14f \t \t \t %3g \n', z');
   else
        fprintf('                 roots ')
        fprintf('   \t\t\t\t\t\t     multiplicities\n');
        fprintf('\n');
        fprintf('%22.15f + %22.15f i \t \t %3g \n', ...
            [real(z(:,1)),imag(z(:,1)),z(:,2)]');
   end;       
