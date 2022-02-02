function [p,z] = jt07d
%
%  test polynomial suggested by Jenkins and Traub
%
   a = 10^(-7);
   p = poly([.001, .01, .1, .1+a*i, .1-a*i, 1, -10]);
   z = [.001, .01, .1, .1+a*i, .1-a*i, 1, -10];
   z = [z',ones(7,1)];
   y = [-10, 1; 1, 1; 0.01, 1; 0.001, 1; 0.1, 3];
   fprintf('\n');
   fprintf(' Ill-conditioned polynomial \n');
   fprintf(' It is constructed using\n')
        fprintf('                 roots ')
        fprintf('   \t\t\t\t\t\t     multiplicities\n');
        fprintf('\n');
        fprintf('%22.15f + %22.15f i \t \t %3g \n', ...
            [real(z(:,1)),imag(z(:,1)),z(:,2)]');
   fprintf('\n');
   fprintf(' However, the polynomial is closer to having \n');
        fprintf('                 roots         multiplicities\n');
        fprintf('\n');
        fprintf('%25.15f \t \t \t %3g \n', y');
