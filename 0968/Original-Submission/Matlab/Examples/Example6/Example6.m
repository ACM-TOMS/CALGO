%
%  Call to disode45
%
 options=disodeset('AbsTol',1.e-4,'RelTol',1.e-4,'ActionSwitch',@actionatswitch6);
 y0=[15;1];
 [tout,yout,tdis,ydis,idis,stats]=disode45(@fun6, @gfun6,[0,20], y0, options);
 plot(tout,yout(:,1));
 xlabel({'t'});
 ylabel({'x(t)'});