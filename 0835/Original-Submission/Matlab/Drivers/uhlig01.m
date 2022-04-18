function [p,z] = uhlig01
%
%  test polynomial used by F. Uhlig
%
   a = 0.01;
   p1 = [1,0,0,0,-a^4];
   p2 = poly([a,a,a,a]);
   p = conv(p1,p2);
   z = [-0.01000000000000 - 0.00000000000000i
 -0.00000000000000 + 0.01000000000000i
 -0.00000000000000 - 0.01000000000000i
  0.01000000000000 - 0.00000000000000i];
   z = [z, [1,1,1,5]'];
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
