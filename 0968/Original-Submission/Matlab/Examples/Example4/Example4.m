%
%  Call to disode45
%
 options=disodeset('Refine',10);
 y0 = [0;0];
 [tout,yout,tdis,ydis,idis,stats]=disode45(@fun4, ...
                                    @gfun4,[0,3], y0,options);
 plot(yout(:,1),yout(:,2),[0.005,0.005],[-0.15,0.25],'r--', ...
                                    [0.005, 0.008],[0,0],'r--');
 xlabel({'y_2(t)'});
 ylabel({'y_1(t)'});