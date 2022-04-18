function scgui(fig)
%SCGUI  Create graphical user interface for SCM Toolbox.
%
%	By itself, SCGUI creates the graphical user interface (GUI)
%	menus for the Schwarz-Christoffel Toolbox in the current
%	figure window.
%
%	SCGUI(FIG) creates the GUI in figure window FIG.
%
%	Use of the GUI is straightforward.  For complete details, see
%	the user's guide.
%
%       See also SCGGET, SCGSET.
%
%	Written by Toby Driscoll.  Last updated 5/31/95.

if nargin < 1
  fig = gcf;
end
figure(fig)
clf
reset(fig)

% menus: column 1 holds menu handles, col 2 is menu group id:
%   0: always available
%   1: available if a polygon has been input
%   2: available if parameter problem has been solved

menus = zeros(30,2);

menus(1) = uimenu('label','Schwarz-Christoffel');

menus(4) = uimenu(menus(1), 'label','Properties...',...
    'call','scguicb(''prop'')','user',[-1,fig,0,8,10,10]);
menus(2) = uimenu(menus(1),'label','Draw new polygon',...
    'separator','on','call','scguicb(''draw'')','interrupt','yes');
menus(5) = uimenu(menus(1),'label','Modify polygon',...
    'call','scguicb(''modify'')','interrupt','yes');
menus(3) = uimenu(menus(1),'label','Load data file...',...
    'call','scguicb(''load'')','interrupt','yes');
menus(5,2) = 1;

menus(10) = uimenu(menus(1),'label','Save data file...',...
    'call','scguicb(''save'')','interrupt','yes');

menus(12) = uimenu(menus(1),'label','Solve parameter problem',...
    'separator','on');
menus(13) = uimenu(menus(12), 'label','half plane -> polygon', ...
    'call','scguicb(''hp2p'')','interrupt','yes');
menus(14) = uimenu(menus(12), 'label','disk -> polygon', ...
    'call','scguicb(''d2p'')','interrupt','yes');
menus(15) = uimenu(menus(12), 'label','disk -> exterior polygon', ...
    'call','scguicb(''d2ep'')','interrupt','yes');
menus(16) = uimenu(menus(12), 'label','strip -> polygon', ...
    'call','scguicb(''st2p'')','interrupt','yes');
menus(17) = uimenu(menus(12), 'label','rectangle -> polygon', ...
    'call','scguicb(''r2p'')','interrupt','yes');
menus(18) = uimenu(menus(12), 'label','continuation', ...
    'call','scguicb(''contin'')','interrupt','yes');
menus([10,12:17],2) = ones(7,1);
menus(18,2) = 2;

menus(20) = uimenu(menus(1),'label','Display results', ...
    'call','scguicb(''disp'')', 'enable','off','interrupt','yes');
menus(21) = uimenu(menus(1),'label','Plot grid image',...
    'call','scguicb(''plot'')', 'enable','off','interrupt','yes');
menus(20:21,2) = 2*ones(2,1);

menus(22) = uimenu(menus(1), 'label','Point source',...
    'call', 'scguicb(''source'')');
menus(22,2) = 2;

% Save menus in figure userdata.
set(fig,'userdata',menus);

% Force compilation of scguicb to speed up later.
if 0, scguicb('draw'), end

scgenable(fig,1:2,'off');


