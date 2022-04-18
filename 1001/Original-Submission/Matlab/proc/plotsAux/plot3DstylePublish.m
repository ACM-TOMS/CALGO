%% plot3DstylePublish
% Define style to publish (seti.plotPublish = 1) figures in case of 3D.
% See also: <setGeomSim.html>.
%
%% Code

% Grid
set(gca, 'XTick', -0.6:0.2:0.6);
set(gca, 'YTick', -0.6:0.2:0.6);
set(gca, 'ZTick', -0.6:0.2:0.6);

% Delete labels
xlabel('');
ylabel('');
zlabel('');

% Grid and box on
grid on;
box on;

% Delete numbers on axis
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
set(gca, 'ZTickLabelMode', 'manual', 'ZTickLabel', []);

% set(gca, 'XTickLabelMode', 'auto', 'XTickLabel','auto');
%set(gca,'FontSize',16);
%xlabel('x_1');
%ylabel('x_2');
%zlabel('x_3');

% Another grid
%set(gca, 'XTick', -0.5:0.25:0.5);
%set(gca, 'YTick', -0.5:0.25:0.5);
%set(gca, 'ZTick', -0.5:0.25:0.5);
