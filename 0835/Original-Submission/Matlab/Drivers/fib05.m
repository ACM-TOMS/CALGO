function [p,z] = fib05
%
%  test polynomial suggested by Goedecker
%
    n = 5;
    p = [-1,ones(1,n)];
    z = [-.6783507129699967-.4585361872731445*i,  1; ...
         -.6783507129699967+.4585361872731445*i,  1; ...
         .19537659464725405-.8488536405462456*i,  1; ...
         .19537659464725405+.8488536405462456*i,  1; ...
         1.9659482366454853,                      1];
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
