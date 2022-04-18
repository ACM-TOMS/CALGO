%
%  Call to disode45
%
 y0=[3;0];
 [tout,yout,tdis,ydis,idis,stats]=disode45(@fun0, @gfun0,[0,20], y0);
 plot(tout,yout(:,1),'k',tout,yout(:,2),'k--',tdis,ydis(:,2),'ro');
 xlabel({'t'});
 ylabel({'y(t)'});
 annotation('textbox',...
    [0.17 0.80 0.0252960172228202 0.0299785867237687],...
    'String','y(t)',...
    'LineStyle','none');
 annotation('textbox',...
    [0.24 0.75 0.0252960172228202 0.0299785867237687],...
    'String','y''(t)',...
    'LineStyle','none');

 figure(2)
 plot(yout(:,1),yout(:,2),'k-',ydis(:,1),ydis(:,2),'ro');
 xlabel({'y(t)'});
 ylabel({'y''(t)'});