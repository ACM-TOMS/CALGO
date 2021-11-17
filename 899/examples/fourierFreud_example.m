
% fourierFreud_example.m 

   clear, home, close all
   N = 384; M = 498;

   h = @(x) 0.*(x<-0.5) + sin(cos(x)).*( (x<=0.5) - (x<-0.5) ) + 0.*(x>0.5);

     x = -1 + 2*(0:N-1)/N;
    xp = -1 + 2*(0:M-1)/M;

    f = h(x);
   fp = h(xp);

    [uc,ak] = fourierInterpolation(f,xp);

    edgeLocations = [-1 -0.5 0.5 1];
    M = [2 12 2];

% automate the selection of m in each subinterval
%       [up,gh,Lambda] = frp(edgeLocations,ak,xp,0,length(ak));    % use spectral coefficients to reconstruct
%       [up,gh,Lambda] = frp(edgeLocations,h,xp,0,length(ak));    % use evenly space function values to reconstruct

%  specify m in each subinterval    
    [up,gh,Lambda] = frp(edgeLocations,ak,xp,0,length(ak),M);    % use spectral coefficients to reconstruct
%     [up,gh,Lambda] = frp(edgeLocations,h,xp,0,length(ak),M);    % use evenly space function values to reconstruct
    
    plot(xp,fp,'r',xp,up,'b')
    xlabel 'x', ylabel 'f(x)'
    axis([-1 1 -0.125 0.9])

    figure
    semilogy(xp,abs( fp - uc ),'g',xp,abs( fp - up ),'r')
    xlabel 'x', ylabel 'log|error|'
