function [p,z] = jt11a
%
%  test polynomial suggested by Jenkins and Traub
%
    m = 15;
    y = exp(i*pi*[(1-m):(m-1)]/(2*m));
    z = 0.9*exp(i*pi*[m:3*m]/(2*m));
    p1 = poly(y); p2 = poly(z);
    p = conv(p1,p2);
    z = [[y';z'],ones(4*m,1)];
    fprintf('\n');
    %fprintf(' Roots are all simple and sensitive, see plot ');
    %plot(real(z(:,1)),imag(z(:,1)),'+');
