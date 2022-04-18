function [p,z] = lgd10
%
%  test polynomial suggested by Goedecker
%   Legendre polynomial of degree 10
%
    p = lgd(10);
    z = [0.97390652851717
   0.86506336668898
   0.67940956829903
  -0.97390652851717
  -0.86506336668899
  -0.67940956829903
   0.43339539412925
  -0.43339539412925
   0.14887433898163
  -0.14887433898163];
    z = [z,ones(10,1)];
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
