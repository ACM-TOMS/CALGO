function [p,z] = jt06
%
%  test polynomial suggested by Jenkins and Traub
%
   p = poly([.1,1.001,.998,1.00002,.99999]);
   z = [.1,1.001,.998,1.00002,.99999];
   z = [z',ones(5,1)];
   fprintf('\n');
   fprintf(' Ill-conditioned polynomial \n');
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
