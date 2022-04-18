%
%  Definition of the vector field
%
 function f=fun1(t,y)
   E=1.0;
   if y(1) > 0
      mu1=1.0;
   else
      mu1=0.6;
   end
   if y(2) > 0
      mu2=0.2;
   else
      mu2=0.5;
   end
   f=[y(3);y(4);-E*(y(1)-y(2))-mu1*sign(y(3)); ...
                    -E*(y(2)-y(1))-mu2*sign(y(4))] ;
 end
