%
%  Definition of the vector field
%
 function f=fun(t,y)
      k=1.0;
      r=0.2;
      a0=1.0;
      w=0.7;
      Fc=0.4;
      v=cos(t)+0.7;
      f=[y(2);-k*y(1)-2*r*y(2)+a0*cos(w*t)-Fc*sign(y(2)-v)];
 end