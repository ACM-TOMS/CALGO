function [p,z] = inex01
%
%  test polynomial suggested by Goedecker
%
    p = poly([(10/11)*[1,1,1,1,1],(20/11)*[1,1,1],...
            (30/11)*[1,1]]);
    p = round(10^7*p)/10^7;
    z = [10/11, 5; 20/11, 3; 30/11, 2];
    z = z(3:-1:1,:);
    fprintf('\n');
    fprintf(' Coefficients are rounded up at the 7-th digits after\n');
    fprintf('    decimal point. Originally, \n');
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
