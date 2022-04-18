
% fourierAndChebyshevGibbs_example.m 

   clear, home, close all
   N = 200; M = 298;

   h = @(x) 0.*(x<-0.5) + sin(cos(x)).*( (x<=0.5) - (x<-0.5) ) + 0.*(x>0.5);

        x = -cos((0:N)*pi/N);
       xp = linspace(-1,1,M);

    f = h(x);
   fp = h(xp);

    [uc,ak] = chebyshevInterpolation(f,xp);

    plot(xp,fp,'r',xp,uc,'g')
    xlabel 'x', ylabel 'f(x)'
    axis([-1 1 -0.125 0.9])
    title 'Chebyshev'
    
    
     x = -1 + 2*(0:N-1)/N;
    xp = -1 + 2*(0:M-1)/M;

    f = h(x);
   fp = h(xp);

    [uc,ak] = fourierInterpolation(f,xp);

    
    figure
    
    plot(xp,fp,'r',xp,uc,'g')
    xlabel 'x', ylabel 'f(x)'
    axis([-1 1 -0.125 0.9])
    title 'Fourier'
    


