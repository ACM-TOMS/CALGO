
% fourierDTV_example.m 

   clear, home, close all
   N = 200; M = 298;

   h = @(x) 0.*(x<-0.5) + sin(cos(x)).*( (x<=0.5) - (x<-0.5) ) + 0.*(x>0.5);

     x = -1 + 2*(0:N-1)/N;
    xp = -1 + 2*(0:M-1)/M;

    f = h(x);
   fp = h(xp);

    [uc,ak] = fourierInterpolation(f,xp);

    timeSteps = 100;
    lambda = 10;
   
    up = digitalTotalVariationFilterPeriodic(uc,lambda,timeSteps); 
     
%    plot(xp,fp,'r',xp,uc,'g',xp,up,'b')
    plot(xp,fp,'r',xp,up,'b')
    xlabel 'x', ylabel 'f(x)'
    axis([-1 1 -0.125 0.9])

    figure
    
    semilogy(xp,abs( fp - uc ),'g',xp,abs( fp - up ),'r')
    xlabel 'x', ylabel 'log|error|'
    axis([-1 1 10e-8 1])
