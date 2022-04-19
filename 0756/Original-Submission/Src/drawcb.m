function [x,y] = drawcb(event,cmd)
%DRAWCB (not intended for calling directly by the user)
%       Callback for DRAWPOLY.
%
%	Written by Toby Driscoll.  Last updated 5/31/95.

global DRP_LINE DRP_PT DRP_AUX

% DRP_LINE: The preview line.  Also, its Userdata stores the points
% selected thus far.

% DRP_PT: Selected point.  The ButtonUpFcn sets this according to what it
% sees in the Userdata of DRP_LINE.  

% DRP_AUX: Parameters.
%   1: Drawing mode.  
%     0=normal, 1=2nd point of Inf vertex, 2=after an Inf vertex
%   2,3: Grid spacing.  If both zero, no grid; else x,y spacings.
%   4: Angle quantization.  If zero, inactive; else fundamental angle.
%   5: Length quantization.  If zero, inactive; else fundamental length.

if strcmp(event,'getpoint')		% get a point from user
  DRP_PT = [];
  set(cmd,'windowbuttonmotionfcn','drawcb(''move'');');
  drawnow
  while isempty(DRP_PT) 		% wait for point selection
    drawnow				
  end;
  x = DRP_PT(1);
  y = DRP_PT(2);

elseif strcmp(event,'move')		% mouse motion
  ptrpos = get(gca,'currentpoint');
  P = ptrpos(1,1:2);
  axlim = axis;
  pts = get(DRP_LINE, 'userdata');
  [m,junk] = size(pts);
  mode = DRP_AUX(1);
  grid = DRP_AUX(2:3)';
  qang = DRP_AUX(4);
  qlen = DRP_AUX(5);

  % Modify point to meet mode constraints.
  if mode==0 				% normal mode
    % No constraints.
  elseif mode==1 			% infinite mode
    % Point may not be inside axes box.
    if all(P>axlim([1,3])) & all(P<axlim([2,4]))
      [junk,j] = min(abs([P(1)-axlim(1:2);P(2)-axlim(3:4)]'));
      P = axlim(j);
    end
  elseif mode==2 			% post-infinite mode
    % Point may not be outside axes box.
    P(1) = min(max(P(1),axlim(1)),axlim(2));
    P(2) = min(max(P(2),axlim(3)),axlim(4));
    % Angle from infinity may not be acute.
    ang = scangle([pts(m-3:m-1,1);P(1)]+i*[pts(m-3:m-1,2);P(2)]);
    ang = ang(2:3);
    ang(ang>0) = ang(ang>0) - 2;
    ang = sum(ang);
    if ang > -1
      % Would be illegal.  Project to make ang=-1.
      A = pts(m-1,:);
      B = pts(m-3,:) + A - pts(m-2,:);
      P = A + ((B-A)*(P-A)')/((B-A)*(B-A)')*(B-A);
      ang = -1;
      qang = 0;				% override other restrictions
      grid = 0;
    elseif ang < -3
      % It's illegal.  Is it even possible?
      P = [NaN,NaN];
      ang = NaN;
    end
  end

  % Modify point to meet angle, length, or grid constraints.

  if any(mode==[0,2]) & qang & (m > 2)	% quantized angle
    % Find arg of new side which meets quantization requirements. 
    if mode==0
      ang = scangle([pts(m-2:m-1,1);P(1)]+i*[pts(m-2:m-1,2);P(2)]);
      ang = qang*round(ang(2)/qang);
      theta = atan2(pts(m-1,2)-pts(m-2,2),pts(m-1,1)-pts(m-2,1))-pi*ang;
    elseif mode==2
      % ang was computed above
      ang = qang*round(ang/qang);
      theta = atan2(pts(m-2,2)-pts(m-3,2),pts(m-2,1)-pts(m-3,1))-pi*ang;
    end
    % Project P to correct angle.
    A = pts(m-1,:);
    BA = [cos(theta),sin(theta)];
    P = A + ((BA)*(P-A)')*(BA);
    grid = 0;
  end
  if (mode==0) & qlen & (m > 1)		% quantized length
    A = pts(m-1,:);
    len = norm(P-A);
    fixlen = qlen*(round(len/qlen));
    P = A + fixlen/len*(P-A);
    grid = 0;
  end
  if any(mode==[0,1,2]) & all(grid)	% snap to grid
    minxy = axlim([1,3]);
    P = minxy + grid.*(round((P-minxy)./grid));
  end
  
  % Update.
  if m > 1
    if ~(mode==1)			% preview line
      set(DRP_LINE, 'xdata',[pts(m-1,1),P(1)], 'ydata',[pts(m-1,2),P(2)]);
    end
    set(DRP_LINE, 'userdata',[pts(1:m-1,:);P]);
  else
    set(DRP_LINE, 'userdata',P);
  end
  drawnow

elseif strcmp(event,'up')		% mouse up
  pts = get(DRP_LINE, 'userdata');
  m = size(pts,1);
  if ~isnan(pts(m,1))			% valid point
    set(gcf,'windowbuttonmotionfcn','');
    set(DRP_LINE,'xdata',[pts(m,1),NaN], 'ydata',[pts(m,2),NaN])
    set(DRP_LINE,'userdata',[pts;[NaN,NaN]])
    DRP_PT = pts(m,:);
  end
  
elseif strcmp(event,'control')		% ui control
  axlim = axis;
  data = get(DRP_LINE,'userdata');
  control = get(gca,'userdata');
  mode = DRP_AUX(1);
  if strcmp(cmd,'g')			% grid feature
    if get(control(2),'value') 		% grid on
      % Turn off quantizations, if now on.
      if get(control(3),'value')
	set(control(3),'value',0)
	drawcb('control','a');
      end
      if get(control(4),'value')
	set(control(4),'value',0)
	drawcb('control','l');
      end
      % Get number of grid points.
      N = round(get(control(5),'value'));
      % Set up allowable x,y values.
      x = linspace(axlim(1),axlim(2),N+1);
      y = linspace(axlim(3),axlim(4),N+1);
      set(gca,'xticklabelmode','auto')
      set(gca,'yticklabelmode','auto')
      set(gca,'xtick',x,'ytick',y)
      % For clarity, keep only about eight of the labels.
      keep = [1,3:ceil(N/8):N-1,N+1];
      xl = get(gca,'xticklabels');
      p = min(size(xl,2),4);
      xlnew = setstr(ones(N+1,1)*blanks(p));
      xlnew(keep,:) = xl(keep,1:p);
      yl = get(gca,'yticklabels');
      p = min(size(yl,2),4);
      ylnew = setstr(ones(N+1,1)*blanks(p));
      ylnew(keep,:) = yl(keep,1:p);
      % Make it so.
      set(gca,'xticklabels',xlnew,'xgrid','on',...
	  'yticklabels',ylnew,'ygrid','on')
      drawnow
      DRP_AUX(2:3) = [x(2)-x(1),y(2)-y(1)];
    else 				% grid off
      set(gca,'xtickmode','auto','xticklabelmode','auto','xgrid','off')
      set(gca,'ytickmode','auto','yticklabelmode','auto','ygrid','off')
      drawnow
      DRP_AUX(2:3) = [0,0];
    end
  elseif strcmp(cmd,'a')		% quantize angle
    DRP_AUX(4) = get(control(3),'value')/round(get(control(6),'value'));
    % Turn off grid, if now on.
    if get(control(2),'value') & get(control(3),'value')
      set(control(2),'value',0)
      drawcb('control','g');
    end
  elseif strcmp(cmd,'l') 		% quantize length
    pct = 1/round(get(control(7),'value'));
    DRP_AUX(5) = get(control(4),'value')*pct*(axlim(2)-axlim(1));
    % Turn off grid, if now on.
    if get(control(2),'value') & get(control(4),'value')
      set(control(2),'value',0)
      drawcb('control','g');
    end
  elseif strcmp(cmd,'sg')
    set(control(8),'string',sprintf('1/%i',round(get(control(5),'value'))));
    drawcb('control','g');
  elseif strcmp(cmd,'sa')
    set(control(9),'string',sprintf('pi/%i',round(get(control(6),'value'))));
    drawcb('control','a');
  elseif strcmp(cmd,'sl')
    set(control(10),'string',sprintf('1/%i',round(get(control(7),'value'))));
    drawcb('control','l');
  end
  % Call move so new restrictions take effect immediately.
  if ~strcmp(cmd(1),'s')
    drawcb('move');
  end

elseif strcmp(event,'key')		% key pressed
  cmd = lower(get(gcf,'currentchar'));
  drawcb('control',cmd);

end



