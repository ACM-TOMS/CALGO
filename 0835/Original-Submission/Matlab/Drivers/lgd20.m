function [p,z] = lgd20
%
%  test polynomial suggested by Goedecker
%   Legendre polynomial of degree 20
%
    p = lgd(20);
    z = [-.9931285991850949, ... 
        -.9639719272779138, ... 
        -.9122344282513259, ... 
        -.8391169718222188, ... 
        -.7463319064601508, ... 
        -.6360536807265150, ... 
        -.5108670019508271, ... 
        -.3737060887154196, ... 
        -.2277858511416451, ... 
        -.7652652113349733e-1, ... 
        .7652652113349733e-1, ... 
        .2277858511416451, ... 
        .3737060887154196, ... 
        .5108670019508271, ... 
        .6360536807265150, ... 
        .7463319064601508, ... 
        .8391169718222188, ... 
        .9122344282513259, ... 
        .9639719272779138, ... 
        .9931285991850949];
    z = [z',ones(20,1)];
    fprintf('Legendre polynomial, roots are sensitive \n');
        fprintf('                 roots         multiplicities\n');
        fprintf('\n');
        fprintf('%25.14f \t \t \t %3g \n', z');
