function [p,z] = farloi01
%
% Farmer-Loizou
%
    p = conv([1,1,2],[1,1,3]);
    p = conv(p,p); p = conv(p,p);
    z = [-0.50000000000000 + 1.65831239517770i
        -0.50000000000000 - 1.65831239517770i
        -0.50000000000000 + 1.32287565553230i
        -0.50000000000000 - 1.32287565553230i];
    z = [z,4*ones(4,1)];
    
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
