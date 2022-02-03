
% edgeDectionChebyshev_example

   clear, home, close all
   N = 200; M = 300;

   h = @(x) cos(0.5*pi.*x).*(x<-0.5) + (x.^3 - sin(1.5*pi.*x) + 1.0).*( (x<=0.5) - (x<-0.5) ) + (x.^2 + 4*x.^3 - 5*x).*(x>0.5);

   x = -cos((0:N)*pi/N);
  xp = linspace(-1,1,M);

   f = h(x);
   fp = h(xp);

   [uc,ak] = chebyshevInterpolation(f,xp);
    plot(xp,fp,'r',xp,uc,'g')

    Q = 2; J = 30; eta=2;


   [S,uE,uN] = edgeDetectChebyshev(ak,J,Q,eta,0);
   ind = find(uN>0); figure
   plot(x,uE,'g',x,uN,'b',x(ind),uN(ind),'k.')
