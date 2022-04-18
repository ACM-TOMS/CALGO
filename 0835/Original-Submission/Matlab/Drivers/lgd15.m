function [p,z] = lgd15
%
%  test polynomial suggested by Goedecker
%   Legendre polynomial of degree 15
%
    p = lgd(15);
    z = [                  0
   0.98799251802063
   0.93727339240040
   0.84820658341064
   0.72441773136009
   0.57097217260857
   0.39415134707756
  -0.98799251802067
  -0.93727339240030
  -0.84820658341075
  -0.72441773136004
   0.20119409399743
  -0.57097217260858
  -0.39415134707756
  -0.20119409399743];
    z = [z,ones(15,1)];
    fprintf('\n');
   if norm(imag(z(:,1))) == 0 
        fprintf('                 roots         multiplicities\n');
        fprintf('\n');
        fprintf('%25.14f \t \t \t %3g \n', z');
   else
        fprintf('                 roots ')
        fprintf('   \t\t\t\t\t\t     multiplicities\n');
        fprintf('\n');
        fprintf('%22.14f + %22.14f i \t \t %3g \n', ...
            [real(z(:,1)),imag(z(:,1)),z(:,2)]');
   end;       
