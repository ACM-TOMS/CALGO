%
%  Definition of the vector field
%
 function f=fun0(t,y)
      Fc=0.4;
      f=[y(2);-y(1)-Fc*sign(y(2))] ;
 end
%
%  Definition of the switching surface
%
 function [g,isterminal,direction]=gfun0(t,y)
      g=y(2);
      isterminal=0;
      direction=0;
 end
 