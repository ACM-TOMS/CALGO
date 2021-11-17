%
%  Definition of the switching surfaces
%
 function [g,isterminal,direction]=gfun2(t,y)
    v=cos(t)+0.7;
    g=y(2)-v;
    isterminal=0;
    direction=0;
 end