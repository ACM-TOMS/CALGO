
%
%  Definition of the switching surfaces
%
 function [g,isterminal,direction]=gfun1(t,y)
   g=[y(3); y(4); y(1); y(2)];
   isterminal=[0;0;0;0];
   direction=[0;0;0;0];
 end