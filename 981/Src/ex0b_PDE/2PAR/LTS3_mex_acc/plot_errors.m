% plot_errors.m
% this script is called by MAIN

if ABSERR
    ERRstr='Abs. Err.';
else
    ERRstr='Rel. Err.';
end

AX = [Xmin Xmax  Tmin Tmax  min([min(min(ERROR1)) min(min(ERROR2)) tol])  max([max(max(ERROR1)) max(max(ERROR2)) tol])];


%% OMP_Talbot11
figure((TOLmenu-1)*4+1); clf; h=surf(n,t,ERROR1); set(h,'EdgeColor',[.5 .5 .5])
    set(gca,'ZScale','log','FontSize',14); xlabel('\eta'); ylabel('t'); zlabel(ZLABEL)
    NV=10; [Vt,Vn]=meshgrid(linspace(Tmin,Tmax,NV),linspace(Xmin,Xmax,NV));
    hold on; h=mesh(Vn,Vt,tol*ones(NV,NV)); set(h,'EdgeColor','k')%'r','LineWidth',2)
    axis(AX); V=axis; colormap('gray'); set(h,'FaceAlpha',0.75)
    if V(end)/tol < 10
        h=text(V(1),V(3),tol,'tol'); set(h,'VerticalAlignment','top','HorizontalAlignment','center','FontWeight','Bold','Color','k') %'r')
    else
        h=text(V(1),V(3),tol,'tol '); set(h,'VerticalAlignment','middle','HorizontalAlignment','right','FontWeight','Bold','Color','k') %'r')
    end
    if V(5) == 0, V(5)=1e-16; axis(V); end
    ZT=[ceil(log10(V(5))) floor(log10(V(6)))]; ZT=unique(sort([round(mean(ZT)) ZT])); ZT=10.^ZT';
    set(gca,'ZTick',ZT,'ZMinorGrid','off','YMinorGrid','off','XMinorGrid','off')
    axis ij; box on
    ht=title({['OMP\_Talbot11\_DE: abs.err. for \eta\in[' num2str(Xmin) ', ' num2str(Xmax) '],  t\in[' num2str(Tmin) ', ' num2str(Tmax) ']'];[repmat(' ',1,45) 'tol=' num2str(tol,'%4.0e') ', NOPTS = ' num2str(NOPTS1)]});
    set(ht,'FontSize',16,'FontWeight','normal')
%{
annotation('textbox',[.02 0 .1 .05],'String',['Example ' EXAMPLE],'FitBoxToText','on','HorizontalAlignment','center')
FIGname=['EX' EXAMPLE '_accuracy_OMP_Talbot11_tol' num2str(TOLmenu)];
set(gcf,'Name',FIGname)
if SAVEflag
    saveas(gcf,[FIGname '.fig']); %saveas(gcf,[FIGname '.png'])
    saveas(gcf,[FIGname '.eps'],'epsc2') % print format = '-depsc2'
end
%}


%% OMP_Talbot12
figure((TOLmenu-1)*4+2); clf; h=surf(n,t,ERROR2); set(h,'EdgeColor',[.5 .5 .5])
    set(gca,'ZScale','log','FontSize',14); xlabel('\eta'); ylabel('t'); zlabel(ZLABEL)
    NV=10; [Vt,Vn]=meshgrid(linspace(Tmin,Tmax,NV),linspace(Xmin,Xmax,NV));
    hold on; h=mesh(Vn,Vt,tol*ones(NV,NV)); set(h,'EdgeColor','k')%'r','LineWidth',2)
    axis(AX); V=axis; colormap('gray'); set(h,'FaceAlpha',0.75)
    if V(end)/tol < 10
        h=text(V(1),V(3),tol,'tol'); set(h,'VerticalAlignment','top','HorizontalAlignment','center','FontWeight','Bold','Color','k') %'r')
    else
        h=text(V(1),V(3),tol,'tol '); set(h,'VerticalAlignment','middle','HorizontalAlignment','right','FontWeight','Bold','Color','k') %'r')
    end
    if V(5) == 0, V(5)=1e-16; axis(V); end
    ZT=[ceil(log10(V(5))) floor(log10(V(6)))]; ZT=unique(sort([round(mean(ZT)) ZT])); ZT=10.^ZT';
    set(gca,'ZTick',ZT,'ZMinorGrid','off','YMinorGrid','off','XMinorGrid','off')
    axis ij; box on
    ht=title({['OMP\_Talbot12\_DE: abs.err. for \eta\in[' num2str(Xmin) ', ' num2str(Xmax) '],  t\in[' num2str(Tmin) ', ' num2str(Tmax) ']'];[repmat(' ',1,45) 'tol=' num2str(tol,'%4.0e') ', NOPTS = ' num2str(NOPTS1)]});
    set(ht,'FontSize',16,'FontWeight','normal')
%{
annotation('textbox',[.02 0 .1 .05],'String',['Example ' EXAMPLE],'FitBoxToText','on','HorizontalAlignment','center')
FIGname=['EX' EXAMPLE '_accuracy_OMP_Talbot12_tol' num2str(TOLmenu)];
set(gcf,'Name',FIGname)
if SAVEflag
    saveas(gcf,[FIGname '.fig']); %saveas(gcf,[FIGname '.png'])
    saveas(gcf,[FIGname '.eps'],'epsc2') % print format = '-depsc2'
end
%}


%% OMP_Talbot13
figure((TOLmenu-1)*4+3); clf; h=surf(n,t,ERROR3); set(h,'EdgeColor',[.5 .5 .5])
    set(gca,'ZScale','log','FontSize',14); xlabel('\eta'); ylabel('t'); zlabel(ZLABEL)
    NV=10; [Vt,Vn]=meshgrid(linspace(Tmin,Tmax,NV),linspace(Xmin,Xmax,NV));
    hold on; h=mesh(Vn,Vt,tol*ones(NV,NV)); set(h,'EdgeColor','k')%'r','LineWidth',2)
    axis(AX); V=axis; colormap('gray'); set(h,'FaceAlpha',0.75)
    if V(end)/tol < 10
        h=text(V(1),V(3),tol,'tol'); set(h,'VerticalAlignment','top','HorizontalAlignment','center','FontWeight','Bold','Color','k') %'r')
    else
        h=text(V(1),V(3),tol,'tol '); set(h,'VerticalAlignment','middle','HorizontalAlignment','right','FontWeight','Bold','Color','k') %'r')
    end
    if V(5) == 0, V(5)=1e-16; axis(V); end
    ZT=[ceil(log10(V(5))) floor(log10(V(6)))]; ZT=unique(sort([round(mean(ZT)) ZT])); ZT=10.^ZT';
    set(gca,'ZTick',ZT,'ZMinorGrid','off','YMinorGrid','off','XMinorGrid','off')
    axis ij; box on
    ht=title({['OMP\_Talbot13\_DE: abs.err. for \eta\in[' num2str(Xmin) ', ' num2str(Xmax) '],  t\in[' num2str(Tmin) ', ' num2str(Tmax) ']'];[repmat(' ',1,45) 'tol=' num2str(tol,'%4.0e') ', NOPTS = ' num2str(NOPTS1)]});
    set(ht,'FontSize',16,'FontWeight','normal')
%{
annotation('textbox',[.02 0 .1 .05],'String',['Example ' EXAMPLE],'FitBoxToText','on','HorizontalAlignment','center')
FIGname=['EX' EXAMPLE '_accuracy_OMP_Talbot13_tol' num2str(TOLmenu)];
set(gcf,'Name',FIGname)
if SAVEflag
    saveas(gcf,[FIGname '.fig']); %saveas(gcf,[FIGname '.png'])
    saveas(gcf,[FIGname '.eps'],'epsc2') % print format = '-depsc2'
end
%}


%{
figure(1); surf(n,t,ft); axis tight; xlabel('\eta'); ylabel('t'); title('true u(\eta,t)')
figure(2); surf(n,t,u1); axis tight; xlabel('\eta'); ylabel('t'); title('classical Talbot: u_1(\eta,t)')
figure(3); surf(n,t,u2); axis tight; xlabel('\eta'); ylabel('t'); title('modified Talbot: u_2(\eta,t)')
%}
