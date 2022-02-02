
% chebyshevPade_example.m 

   clear, home, close all
   N = 200; M = 298;

   h = @(x) 0.*(x<-0.5) + sin(cos(x)).*( (x<=0.5) - (x<-0.5) ) + 0.*(x>0.5);

        x = -cos((0:N)*pi/N);
       xp = linspace(-1,1,M);

    f = h(x);
   fp = h(xp);

    [uc,ak] = chebyshevInterpolation(f,xp);

  
    
    M = 14;
    Nc = 0;
    up = chebyshevPade(ak,M,Nc,xp)';
    
    
%    plot(xp,fp,'r',xp,uc,'g',xp,up,'b')
    plot(xp,fp,'r',xp,up,'b')
    xlabel 'x', ylabel 'f(x)'
    axis([-1 1 -0.125 0.9])

    figure
    semilogy(xp,abs( fp - uc ),'g',xp,abs( fp - up ),'r')
    xlabel 'x', ylabel 'log|error|'
