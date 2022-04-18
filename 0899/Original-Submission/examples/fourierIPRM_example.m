
% fourierIPRM_example.m 

   clear, home, close all
   N = 200; M = 299;

   h = @(x) 0.*(x<-0.5) + sin(cos(x)).*( (x<=0.5) - (x<-0.5) ) + 0.*(x>0.5);

     x = -1 + 2*(0:N-1)/N;
    xp = -1 + 2*(0:M-1)/M;

    f = h(x);
   fp = h(xp);

    [uc,ak] = fourierInterpolation(f,xp);

    edgeLocations = [-1 -0.5 0.5 1];
    M = [2 6 2];
    
    [up,condW] = inverseReprojection(edgeLocations,M,xp,h);     % use evenly space function values to reconstruct
%    [up,condW] = inverseReprojection(edgeLocations,M,xp,ak);    % use spectral coefficients to reconstruct
     
    plot(xp,fp,'r',xp,up,'b')
    xlabel 'x', ylabel 'f(x)'
    axis([-1 1 -0.125 0.9])

    figure
    semilogy(xp,abs( fp - uc ),'g',xp,abs( fp - up ),'r')
    xlabel 'x', ylabel 'log|error|'
    axis([-1 1 10e-8 1])
