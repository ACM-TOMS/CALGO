function [p,z] = fibsq08
%
%  test polynomial suggested by Goedecker
%
    n = 8;
    p = [-1,ones(1,n)];
    p = conv(p,p);
    z = [-.8762862300182460, ... 
        -.6416053894725764-.6063952206057722*i, ... 
        -.6416053894725764+.6063952206057722*i, ... 
        -.4694032461065018e-1-.9030234661229178*i, ... 
        -.4694032461065018e-1+.9030234661229178*i, ... 
        .6286732392246423-.7084725692273577*i, ... 
        .6286732392246423+.7084725692273577*i, ... 
        1.996031179735415];
    z = [z',2*ones(n,1)];
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
