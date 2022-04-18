%% plotExpSetup
% Plot and save experimental set-up in figure 1.
% (i.e. positions of transmitters and receivers).
%
%% Syntax
%
%   plotExpSetup(seti,out)
%
%% Input Arguments
%
% * seti.incPnts  : positions of transmitters
% * seti.measPnts : positions of receivers
% * out           : output depth: no figure (0), plot figure (1), 
%                   plot and save figure (2).
%
% Further details of the fields in the structural array |seti| are given in <expSetup.html>.
%
% Specific fields in |seti| influencing the figures are described in 
% <setGeomSim.html> in the section "Subfunction: setFigureSettings".
%
%% Output Arguments
%
% plot 1, see <start.html>.
%
%% See Also
%
% * <start.html>
% * <expSetup.html>
% * <setGeomSim.html>
% * <savePngFig.html>
%
%% Code

function plotExpSetup(seti,out)

%% 
% *figure 1: source and measurement points plot...*

if out >= 1
    
    figure(1);
    set(gcf,'Visible',seti.plotVisible);
    IP = seti.incPnts; % temporary for short notation
    MP = seti.measPnts;
    if seti.dim == 2
        [~,rInc]  = cart2pol(IP(1,:),IP(2,:));
        [~,rMeas] = cart2pol(MP(1,:),MP(2,:));
    elseif seti.dim == 3
        %
        [~,rInc,~]  = cart2pol(IP(1,:),IP(2,:),IP(3,:));
        [~,rMeas,zMeas] = cart2pol(MP(1,:),MP(2,:),MP(3,:));
        rMeas = max(rMeas,zMeas);
    end
   
    if seti.dim == 2
        % style start
        set(gca,'LooseInset',get(gca,'TightInset')) % weniger Rand
        incMarkerSize = 30;
        measMarkerSize = 7;
        % measLineWidth = 3;
        clf;
        hold on;
        % 2D case
        plot(IP(1,:),IP(2,:),'b.','MarkerSize',incMarkerSize);
        % plot(MP(1,:),MP(2,:),'ro','MarkerSize',measMarkerSize,'LineWidth',measLineWidth);
        plot(MP(1,:),MP(2,:),'rs','MarkerSize',measMarkerSize,'MarkerFaceColor','red');
        clear IP MP;
        factor = max(max(rInc),max(rMeas))*1.2;
        axis(factor*[-1 1 -1 1]); axis square;
        axis off; box off;
        hold off;
        % style end
    elseif seti.dim == 3
        % 3D case
        s = 100; % size
        scatter3(IP(1,:),IP(2,:),IP(3,:),s,'filled','b'); hold on;
        scatter3(MP(1,:),MP(2,:),MP(3,:),s,'filled','rs'); hold off;
        set(gca,'FontSize',seti.pubFontSize);
        % Delete labels
%        xlabel(''); ylabel(''); zlabel('');
        % Grid and box on
        grid on; box on;
        % Delete numbers on axis
%         set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
%         set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
%         set(gca, 'ZTickLabelMode', 'manual', 'ZTickLabel', []);
    end
    if seti.plotPublish
    else
        title('transmitters (blue) and receivers (red)');
    end
    axis square;
    
    if out == 2
        savePngFig(1,0,seti); % set iOut = 0
    end
end

end
