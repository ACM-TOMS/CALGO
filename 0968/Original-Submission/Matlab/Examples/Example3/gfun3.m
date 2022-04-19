%
%  Definition of the switching surfaces
%
 function [g,isterminal,direction]=gfun3(t,y)
    m1=2;
    m2=1;
    r1=0.6;
    r2=0.6;
    k1=30;
    k2=20;
    a0=30;
    b0=35;
    d=1.0;
    ee=0.7;
    w=1.38;
    u1=-(r1/m1)*y(3)-(k1/m1)*y(1)+b0/m1+(a0/m1)*cos(w*t);
    u2=-(r2/m2)*y(4)-(k2/m2)*y(2);
    g=[y(1)-y(2)-d/2;y(2)-y(1)-d/2;u1-u2];
    isterminal=[-1;-1;0];
    direction=[1;1;0];
 end