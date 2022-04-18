function [w,beta,handles] = drawpoly(fig,axlim)
%DRAWPOLY Draw a polygon with the mouse.
%       [W,BETA] = DRAWPOLY allows the user to draw a polygon with the
%       mouse.  Use the mouse to position the crosshair and press the
%       left mouse button to create a vertex.  For use with other S-C
%       Toolbox functions, the vertices must be specified in a
%       "positively oriented" manner; i.e.  counterclockwise for
%       interior polygons and clockwise for exterior regions.  There are
%       several GUI elements added to the figure to help you snap
%       vertices to a grid, get specfic angles, etc.  For the last
%       vertex, use the middle or right mouse button, or double click.
%       Upon return, W is a vector of complex vertices and BETA is a
%       vector of turning angles.
%
%       [W,BETA] = DRAWPOLY(FIG) draws in figure FIG.  [W,BETA] =
%       DRAWPOLY(FIG,AXLIM) also uses AXLIM for the axes limits.
%
%       [W,BETA,H] = DRAWPOLY also returns a vector of handles to the
%       plotted edges.
%
%	See the user's guide for full details.
%	
%	See also PLOTPOLY, MODPOLY, SCGUI.
%
%	Written by Toby Driscoll.  Last updated 5/24/95.

if nargin < 2
  axlim = [-4,4,-4,4];
  if nargin < 1
    fig = gcf;
  end
end

figure(fig);

% See DRAWCB for documentation of globals.
global DRP_LINE DRP_PT DRP_AUX

% Set up figure, axes, etc.
if ~ishold, cla, end
view(2)
oldbuf = get(fig,'windowbuttonupfcn');
oldbdf = get(fig,'windowbuttondownfcn');
set(gca,'xlim',axlim(1:2),'ylim',axlim(3:4),'aspect',[1,NaN],'box','on')
hold on
ptr = get(fig, 'pointer');
set(fig, 'pointer','crosshair');
DRP_LINE = line(NaN,NaN,'linestyle','--','erasemode','xor',...
    'clipping','off');

% Set up check boxes.
oldun = get(fig,'units');
set(fig,'units','centi');
figpos = get(fig,'position');
set(gca,'units','centi')
axpos = get(gca,'pos');
set(gca,'pos',[axpos(1),axpos(2)+1.8,axpos(3:4)])
set(fig,'pos',[figpos(1),figpos(2)-1.8,figpos(3),figpos(4)+1.8])
control(1) = uicontrol('style','frame','units','centi',...
    'pos',[0,0,figpos(3),1.8]);
set(control(1),'units','norm')
offset = max(0,(figpos(3)-11.6)/2);
control(2) = uicontrol('style','check','string','Snap to grid ',...
    'units','centi','pos',[offset .8 3.7 .8]);
set(control(2),'units','norm','call','drawcb(''control'',''g'');')
control(3) = uicontrol('style','check','string','Quantize angle',...
    'units','centi','pos',[offset+3.95 .8 3.7 .8]);
set(control(3),'units','norm','call','drawcb(''control'',''a'');')
control(4) = uicontrol('style','check','string','Quantize length',...
    'units','centi','pos',[offset+7.9 .8 3.7 .8]);
set(control(4),'units','norm','call','drawcb(''control'',''l'');')
% Create sliders.
control(5) = uicontrol('style','slider','min',4,'max',32,...
    'value',16,'units','cent','pos',[offset+.45 .2 1.75 .4]);
set(control(5),'units','norm','call','drawcb(''control'',''sg'');')
control(6) = uicontrol('style','slider','min',2,'max',24,...
    'value',12,'units','cent','pos',[offset+4.4 .2 1.75 .4]);
set(control(6),'units','norm','call','drawcb(''control'',''sa'');')
control(7) = uicontrol('style','slider','min',3,'max',20,...
    'value',8,'units','cent','pos',[offset+8.35 .2 1.75 .4]);
set(control(7),'units','norm','call','drawcb(''control'',''sl'');')
% Text to accompany sliders.
control(8) = uicontrol('style','text','string','1/16',...
    'units','cent','pos',[offset+2.5 .2 1 .4]);
set(control(8),'units','norm')
control(9) = uicontrol('style','text','string','pi/12',...
    'units','cent','pos',[offset+6.3 .2 1 .4]);
set(control(9),'units','norm')
control(10) = uicontrol('style','text','string','1/8',...
    'units','cent','pos',[offset+10.15 .2 1 .4]);
set(control(10),'units','norm')
set(gca,'userdata',control)
drawnow

% Preparation.
DRP_AUX = zeros(5,1);
set(DRP_LINE, 'userdata',[NaN,NaN]);
if ~strcmp(computer,'SUN4')
  % Kludge.  Draw preview line when button is pressed.
  set(fig,'windowbuttondownfcn','drawcb(''move'');');
end
set(fig,'windowbuttonupfcn', 'drawcb(''up'');');
%%set(fig,'keypressfcn','drawcb(''key'');');
drawnow

% Get first vertex.
[x,y] = drawcb('getpoint',fig);
w = x+i*y;
vertices = plot(x,y,'.', 'markersize',12); 
set(DRP_LINE,'xdata',[x NaN], 'ydata',[y,NaN]);
drawnow
button = 1;
n = 1;

% Get rest of vertices.
DRP_PT = [];				% no point selected
mode = 0;				% "normal" mode
while button==1				% until last was selected
  DRP_AUX(1) = mode;
  n = n + 1;
  [x0,y0] = drawcb('getpoint',fig);
  m = length(x);
  x = [x;x0];  y = [y;y0];
  edges(n-1) = plot(x(m:m+1),y(m:m+1),'-');
  set(vertices, 'xdata',x, 'ydata',y);
  drawnow
  if x0>=axlim(1) & x0<=axlim(2) & y0>=axlim(3) & y0<=axlim(4)
    % Finite vertex.
    w(n) = x0+i*y0;
    mode = 0;
  else					% infinite vertex
    % Get re-entry point for next edge.
    DRP_AUX(1) = 1;			% "inf" mode
    set(fig,'pointer','cross')
    [x0,y0] = drawcb('getpoint',fig);
    set(fig,'pointer','crosshair')
    x = [x;x0];  y = [y;y0];
    w(n) = Inf;
    mode = 2;				% "post-inf" mode
  end  % if vertex is infinite
  if n > 2 				% angle at previous vertex
    if ~isinf(w(n-1))
      ang = scangle(x(m-1:m+1)+i*y(m-1:m+1));
      beta(n-1) = ang(2);
    else
      ang = scangle(x(m-2:m+1)+i*y(m-2:m+1));
      ang = ang(2:3);
      ang(ang>0) = ang(ang>0) - 2;
      beta(n-1) = sum(ang);
    end
  end  % if n > 2
  % What kind of button press?
  if strcmp(get(fig,'selectiontype'),'normal')
    button = 1;
  else
    button = 2;
  end
end  % while button==1

m = length(x);
edges(n) = plot(x([m,1]),y([m,1]),'-');
drawnow
% Angle at vertex n.
if ~isinf(w(n))
  ang = scangle(x([m-1:m,1])+i*y([m-1:m,1]));
  beta(n) = ang(2);
else
  ang = scangle(x([m-2:m,1])+i*y([m-2:m,1]));
  ang = ang(2:3);
  ang(ang>0) = ang(ang>0) - 2;
  beta(n) = sum(ang);
end

% Angle at first vertex (necessarily finite).
ang = scangle(x([m,1:2])+i*y([m,1:2]));
beta(1) = ang(2);

% Prepare outputs.
w = w(:);
beta = beta(:);

% Clean up the mess.
set(fig, 'pointer',ptr, 'windowbuttonupfcn',oldbuf,...
    'windowbuttondownfcn',oldbdf,...
    'windowbuttonmotionfcn', '',...
    'pos',figpos,'units',oldun);
set(gca,'pos',axpos,'units','norm')
for j = 1:length(control)
  delete(control(j))
end
delete(DRP_LINE)
clear DRP_LINE DRP_PT DRP_AUX

hold off  
axis auto
handles = plotpoly(w,beta);

