%
%  Call to disode45
%
 options=disodeset('RelTol',1.e-5,'AbsTol',1.e-5, ...
                        'ActionSwitch',@actionatswitch3);
 y0 = [0.2;0.3;0;0];
 [tout,yout,tdis,ydis,idis,stats]=disode45(@fun3, @gfun3, ...
                                       [0,10], y0,options);
 plot(tout,yout(:,1),'k',tout,yout(:,2),'b',tdis(idis<0), ...
         ydis(idis<0,1),'ro',tdis(idis<0),ydis(idis<0,2),'ro')
 xlabel({'t'});
 ylabel({'y_1(t), y_2(t)'});
 annotation('textbox',...
    [0.2 0.86 0.0252960172228202 0.0299785867237687],...
    'String','y_1(t)',...
    'LineStyle','none');
 annotation('textbox',...
    [0.27 0.22 0.0252960172228202 0.0299785867237687],...
    'String','y_2(t)',...
    'LineStyle','none');
 figure(2);
 plot(tout,yout(:,1)-yout(:,2),'k')
 xlabel({'t'});
 ylabel({'y_1(t)-y_2(t)-0.5'});
