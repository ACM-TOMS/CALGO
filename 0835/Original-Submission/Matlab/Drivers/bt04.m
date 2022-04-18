function [p,z] = bt04
%
% Brugnano and Trigiante
%
    p = poly([1,1,1, -1*[1,1,1,1], (.5+i)*[1,1,1], (.5-i)*[1,1,1], ...
            .5*(1+i)*[1,1], .5*(1-i)*[1,1]]);
    z = [1,3; -1,4; .5+i,3; .5-i,3; .5*(1+i),2; .5*(1-i),2];
    
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
