%
%  Call to disode45
%
 options=disodeset('AbsTol',1.e-4,'RelTol',1.e-4, ...
                 'Refine',10, 'ActionSwitch',@actionatswitch5);
 y0=[10; 0];
 [tout,yout,tdis,ydis,idis,stats]=disode45(@fun5, ...
                              @gfun5,[0,20], y0, options);
 plot(tout,yout(:,1));
 xlabel({'t'});
 ylabel({'x(t)'});