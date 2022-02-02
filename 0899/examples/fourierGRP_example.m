
% fourierGRP_example.m 

   clear, home, close all
   N = 200; M = 298;

   h = @(x) 0.*(x<-0.5) + sin(cos(x)).*( (x<=0.5) - (x<-0.5) ) + 0.*(x>0.5);

     x = -1 + 2*(0:N-1)/N;
    xp = -1 + 2*(0:M-1)/M;

    f = h(x);
   fp = h(xp);

    [uc,ak] = fourierInterpolation(f,xp);
    
%    edgeLocations = [-1 -0.48 0.5 1];        % bad example, edge not correctly located
%    lambda = [4 7 4];
%    M = [13 13 13];

    edgeLocations = [-1 -0.5 0.5 1];   % good example
    lambda = [20 25 20];
    M = [2 12 2];
    
    [up,gh] = grp(edgeLocations,lambda,M,ak,xp,0);               % use spectral coefficients to reconstruct
%    [up,gh] = grp(edgeLocations,lambda,M,h,xp,0);             % use evenly space function values to reconstruct
    
    
%    plot(xp,fp,'r',xp,uc,'g',xp,up,'b')
    plot(xp,fp,'r',xp,up,'b')
    xlabel 'x', ylabel 'f(x)'
    axis([-1 1 -0.125 0.9])

    figure
    semilogy(xp,abs( fp - uc ),'g',xp,abs( fp - up ),'r')
    xlabel 'x', ylabel 'log|error|'
