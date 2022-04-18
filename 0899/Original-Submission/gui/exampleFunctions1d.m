% EXAMPLEFUNCTIONS1D  example functions available from the Examples1d menu of the MPT GUI
%
% Inputs:
%   x[]      the grid the function is know on
%  xp[]      the grid the function is interpolated to (identical to x in the PDE case)
%   fCh      function number
% Outputs:
%    f[]    f(x)
%   fp[]    f(xp)
%      h    a handle to the selected function
% Called by:
%    1) exampleFunctionSetup1d.m
% Last modified: October 17, 2007

function [f,fp,h] = exampleFunctions1d(x,xp,fCh)

   switch fCh
     case 0              %  non-periodic
       h = @(x) x;
     case 1              %  step function, non-periodic, jump at x=0
       h = @(x) -1.*(x<=0) + 1.*(x>0);
     case 2               % 
       h = @(x) cos(0.5*pi.*x).*(x<-0.5) + (x.^3 - sin(1.5*pi.*x) + 1.0).*( (x<=0.5) - (x<=-0.5) ) + (x.^2 + 4*x.^3 - 5*x).*(x>0.5);
     case 3
       h = @(x) sin(0.5*pi.*(x+1)).*(x<=0) - sin(0.5*pi.*(x+1)).*(x>0);
     case 4
       h = @(x) 1.*(x>= -0.4 & x<=-0.2) + 1.*(x>= 0.2 & x<=0.4) + ((15/4)*x + 1).*(x>= -0.2 & x<0) + ((-15/4)*x + 1).*(x>= 0 & x<=0.2);
     case 5               % Sharp peak, center
       h = @(x) (-1-x).*(x<=0) + ((1-x).^6).*(x>0);
     case 6               % center step
       K = 0.4;
       h = @(x) 0.*(x<-K) + 1.*( (x<=K) - (x<-K) ) + 0.*(x>K);
     case 7
       h = @(x) exp( cos(8*x.^3 + 1) );
     case 8
       a = 0.5; z = -0.7; delta = 0.005; alpha = 10;  beta = log(2)/(36*delta^2);
       h = @(x) (1/6)*( G(x,beta,z-delta) + G(x,beta,z+delta) + 4*G(x,beta,z) ).*(x>=-0.8 & x<=-0.6) + ...
           1.*(x>=-0.4 & x<=-0.2) + (1 - abs(10*(x-0.1))).*(x>=0 & x<=0.2) + ...
          (1/6)*( F(x,alpha,a-delta) + F(x,alpha,a+delta) + 4*F(x,alpha,a) ).*(x>=0.4 & x<=0.6);
     case 9              % sawtooth
       h = @(x) (x+1.0).*(x<0) + (x-1.0).*(x>=0);
     case 10             % off center sharp peak
       h = @(x) ((2*exp(2*pi.*(x+1))-1-exp(pi))./(exp(pi)-1)).*(x<-0.5) + (-sin(2*pi*x/3 + pi/3)).*(x>=-0.5);
     case 11             %  Shepp-Logan slice
      h = @(x) shepLogan_mpt(x,0.5);
% h = @(x) shepLogan_mpt(x,-0.59);
     case 12             % sin(cos)*I[-0.5,0.5]
       h = @(x) 0.*(x<-0.5) + sin(cos(x)).*( (x<=0.5) - (x<-0.5) ) + 0.*(x>0.5);
     case 13             % discontinuous 1st derivative, non-periodic
       h = @(x) (-x).*(x<0) + (2*x).*(x>=0);
     case 14             % discontinuous 1st derivative, periodic
       h = @(x) abs(x);
     case 15             % analytic and periodic
       h = @(x) exp( sin(3*x*pi) + cos(x*pi) );
     case 16             % discontinuous f(x) and f'(x)
       h = @(x) (-x).*(x<0) + (2*x).*(x>=0 & x<0.5) + 0.*(x>=0.5);
   end

    f = h(x);
   fp = h(xp);


% -------------------- nested functions ---------------------------------


   function r = G(x,beta,z)
      r = exp(-beta.*(x-z).^2);
   end

   function r = F(x,alpha,a)
       r = sqrt( max(1-alpha^2.*(x-a).^2,0) );
   end

% -----------------------------------------------------------------------

end
