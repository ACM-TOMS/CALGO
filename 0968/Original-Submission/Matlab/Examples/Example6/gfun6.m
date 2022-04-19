%
%  Definition of the switching surfaces
%
 function [g,isterminal,direction]=gfun6(t,y)
      g=[y(1)-23.5; y(1)-22];
      isterminal=[-1;-1];  %  Call to gwhenswitch when found
      direction=[1;-1]; % From negative to positive the first one
 end