%
%  Call to disode45
%
 y0 = [3;0];
 [tout,yout,tdis,ydis,idis,stats]=disode45(@fun2, @gfun2,[0,30], y0);
 plot(tout,yout(:,1),'k',tout,yout(:,2),'k--', tdis, ydis(:,1),'ro');
 xlabel({'t'});
 ylabel({'y(t), y''(t)'});
 annotation('textbox',...
    [0.15 0.86 0.0252960172228202 0.0299785867237687],...
    'String','y(t)',...
    'LineStyle','none');
 annotation('textbox',...
    [0.245 0.60 0.0252960172228202 0.0299785867237687],...
    'String','y''(t)',...
    'LineStyle','none');
 figure(2);
 plot(tout,yout(:,2)-0.7-cos(tout)',[0,30],[0,0],'r');
 xlabel({'t'});
 ylabel({'y''(t)-v(t)'});
