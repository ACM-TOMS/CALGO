
% spectralFilter_example

   clear, home, close all
   N = 200; M = 298;

   h = @(x) 0.*(x<-0.5) + sin(cos(x)).*( (x<=0.5) - (x<-0.5) ) + 0.*(x>0.5);

     x = -1 + 2*(0:N-1)/N;
    xp = -1 + 2*(0:M-1)/M;

    f = h(x);
   fp = h(xp);

    [uc,ak] = fourierInterpolation(f,xp);

    filterChoice = 3;
    filterOrder = 4;

    up = filterFourier(ak,xp,filterChoice,filterOrder);
    plot(xp,fp,'r',xp,up,'b')
%    plot(xp,fp,'r',xp,uc,'g',xp,up,'b')
    xlabel 'x', ylabel 'f(x)'
    axis([-1 1 -0.125 0.9])

    figure
    semilogy(xp,abs( fp - uc ),'g',xp,abs( fp - up ),'r')
    axis([-1 1 10e-10 1])
    xlabel 'x', ylabel 'log|error|'
