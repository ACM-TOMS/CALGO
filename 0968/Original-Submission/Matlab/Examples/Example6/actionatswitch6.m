%
%  Definition of the action at switch function
%
 function ysw=actionatswitch6(t,y)
    if y(2)==1,
       ysw=[y(1); -1];
    else
       ysw=[y(1); 1];
    end
 end