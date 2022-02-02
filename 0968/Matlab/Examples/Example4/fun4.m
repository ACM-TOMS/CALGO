function f=fun4(t,y)
      k=210.125;
      c=2.47e+6;
      nu=0.005;
      r=2*sin(14*t);
      if y(1)>nu
        if y(2)>0
           u=c*(y(1)-nu)^(3/2) + 1.98*sqrt(2*c*sqrt(y(1)-nu))*y(2);
        else
           u=c*(y(1)-nu)^(3/2);
        end
      else
        u=0;
      end
      f=[y(2); (-4.1*y(2)-k*y(1)-u-r)/2];
 end