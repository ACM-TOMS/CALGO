%
%  Call to disode45
%
 y0 = [-2;3;0;0];
 [tout,yout,tdis,ydis,idis,stats]=disode45(@fun1, @gfun1,[0,20], y0);
 plot(tout,yout(:,1),'k',tout,yout(:,2),'b', tout, ...
   yout(:,3),'k--',tout,yout(:,4),'b--',tdis,ydis(:,2),'ro');
 xlabel({'t'});
 ylabel({'y_1(t)'});
 annotation('textbox',...
    [0.13 0.86 0.0252960172228202 0.0299785867237687],...
    'String','y_2(t)',...
    'LineStyle','none','Color','b');
 annotation('textbox',...
    [0.245 0.75 0.0252960172228202 0.0299785867237687],...
    'String','y_1(t)',...
    'LineStyle','none');