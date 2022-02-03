 function [g,isterminal,direction]=gfun4(t,y)
     if y(1)<0.005
          g=[y(1)-0.005;1]; 
     else
         g=[y(1)-0.005;y(2)]; 
     end
    isterminal=[0;0];
    direction=[0;-1];
 end