function [p,z] = miyak08
%
%  test polynomial suggested by Goedecker
%  square of Fibocacci polynomial
%
    p = poly([(1.1+1.1*i)*ones(1,4),(3.2+2.3*i)*[1,1],2.1+1.5*i]);
    p = poly([(1.1+1.1*i)*ones(1,32),(3.2+2.3*i)*ones(1,16), ...
            (2.1+1.5*i)*ones(1,8)]);
    z = [[(1.1+1.1*i),(3.2+2.3*i),2.1+1.5*i].',8*[4,2,1]'];
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
