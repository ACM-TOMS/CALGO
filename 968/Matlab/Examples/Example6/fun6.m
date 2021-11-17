%
%  Definition of the vector field
%
 function ydot=fun6(t,y)
     if y(2)==-1
       ydot=[-0.1*(y(1)-18); 0] ;
     else
       ydot=[-0.1*(y(1)-18)+2; 0] ;
     end
 end