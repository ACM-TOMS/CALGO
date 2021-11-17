%
%  Definition of the action at switch function
%
  function ysw=actionatswitch3(t,y)
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
    ysw=y;
    if abs(y(3)-y(4))>1.e-6,
       ysw(3)=((m1-m2*ee)/(m1+m2))*y(3)+((1+ee)*m2/(m1+m2))*y(4);
       ysw(4)=((m2-m1*ee)/(m1+m2))*y(4)+((1+ee)*m1/(m1+m2))*y(3);
       if y(1)-y(2)-d/2>=0
          ysw(1)=y(2)+d/2;
       elseif y(2)-y(1)>=d/2,
          ysw(2)=y(1)+d/2;
       end
    else
       ysw(3)=ysw(4);
       if y(1)-y(2)>=d/2
          ysw(1)=y(2)+d/2;
       elseif y(2)-y(1)>=d/2,
          ysw(2)=y(1)+d/2;
       end
    end