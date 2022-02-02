% EXAMPLEFUNCTIONS2D  example functions available from the Examples2d menu of the MPT GUI
%
% Inputs:
%   (X[][],Y[][])      the grid the function is know on
%   (XP[][],YP[][])    the grid the function is interpolated to (identical to x in the PDE case)
%   fCh                function number
% Outputs:
%    f[][]    f(X,Y)
%   fp[][]    f(XP,YP)
% Called by:
%    1) exampleFunctionSetup2d.m
% Last modified: October 17, 2007

function [f,fp] = exampleFunction2d(X,Y,XP,YP,fCh) 

   switch fCh
     case 1            %  square function          
       h = @(X,Y) 0.*(abs(X)>0.5 & abs(Y)>0.5 ) + 1.*(abs(X)<=0.5 & abs(Y)<=0.5 );
     case 2            % circle non-constant, compact support      
       h = @(X,Y) ((X.^2 + Y.^2)<0.25).*(3*X + 2*Y.^2 + 3);
     case 3            % modified Shepp-Logan brain phantom        
       h = @(X,Y) shepLogan_mpt(X,Y);
     case 4            %  circle, constant, compact support
       h = @(X,Y) 0.*((X.^2 + Y.^2)>0.25) + 1.*((X.^2 + Y.^2)<=0.25);
     case 5            % circle, non-constant, non-compact support 
       h = @(X,Y) ( 10*X + 5 + X.*Y + cos(2*pi.*X.^2) - sin(2*pi.*X.^2) ).*((X.^2 + Y.^2)<=0.25) + ( 10*X - 5 + X.*Y + cos(2*pi.*X.^2) - sin(2*pi.*X.^2) ).*((X.^2 + Y.^2)>0.25);
     case 6
       h = @(X,Y) ((abs(X) + abs(Y))<=0.5).*sin(pi*(X + Y)) - ((abs(X) + abs(Y))>0.75).*sin(pi*(X + Y));
     case 7
       h = @(X,Y)  ( 3.0).*( X<0 & Y<0 & (abs(X) + abs(Y))>=1) + ( 2.0).*( X<0 & Y<0 & (abs(X) + abs(Y))<1) + ...
                   + 2.*( X<0 & Y>=0  ) + 2.*( X>=0 & Y<0  ) + 1.0.*( X>0 & Y>0 & (X.^2 + Y.^2)<0.5  );
   end
   
    f = h(X,Y);
   fp = h(XP,YP);
