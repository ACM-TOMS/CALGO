function scguicb(func)
%SCGUICB (not intended for calling directly by the user)
%       Callback function for SCGUI.
%
%	Written by Toby Driscoll.  Last updated 5/31/95.

fig = gcf;
maptypes = str2mat('hp2p','d2p','d2ep','st2p','r2p');
prefixes = str2mat('hp','d','de','st','r');

if strcmp(func,'draw')			%** draw new polygon
  scgset(fig, 'clear');
  % If drawpoly dies, clear the figure's callbacks.
  errfun = ['set(fig,''windowbuttondownfcn'',''''),',...
          'set(fig,''windowbuttonmotionfcn'',''''),',...
          'set(fig,''windowbuttonupfcn'',''''),',...
          'set(fig,''keypressfcn'',''''),',...
          'return'];
  w = [];
  beta = [];
  eval('[w,beta] = drawpoly(fig);',errfun);
  scgset(fig,'vertices',w,'angles',beta)
  scgenable(fig,2,'off');
  scgenable(fig,1,'on');
  
elseif strcmp(func,'modify')		%** modify polygon
  [w,beta] = scgget(fig,'ver','ang');
  [w,beta] = modpoly(w,beta);
  scgset(fig,'ver',w,'ang',beta)
  scgenable(fig,2,'off')

elseif strcmp(func,'load') 		%** load data file
  [fname,pname] = uigetfile('*.mat','Load S-C data file');
  if fname
    w = []; beta = []; z = []; c = []; maptype = [];
    eval(['load ',pname,fname]);
    scgset(fig, 'clear');
    scgset(fig,'ver',w,'ang',beta,'pre',z,'con',c,'map',maptype);
    if ~isempty(w)
      scgenable(fig,1,'on');
      hold off
      plotpoly(w,beta);
      if ~isempty(z)
	scgenable(fig,2,'on');
      else 
	scgenable(fig,2,'off');
      end
    else 
      scgenable(fig,1,'off');
    end
    if strcmp(maptype,'r2p') 		% rectangle is special
      menus = get(fig,'userdata');
      set(menus(17),'userdata',L);
    end
  end

elseif strcmp(func,'save')		%** save data file
  [fname,pname] = uiputfile('*.mat','Save S-C data file');
  if fname
    [w,beta,z,c,maptype] = scgget(fig); 
    evalstr = ['save ',pname,fname,' w beta z c maptype'];
    if strcmp(maptype,'r2p') 		% rectangle is special
      menus = get(fig,'userdata');
      L= get(menus(17),'userdata');
      evalstr = [evalstr,' L'];
    end
    eval(evalstr);
  end

elseif strcmp(func,'hp2p') 		%** half plane -> polygon
  [w,beta] = scgget(fig, 'vertices','angles');
  [w,beta] = scfix('hp',w,beta);
  scgset(fig,'ver',w,'ang',beta)
  disp('Solving parameter problem...')
  [trace,tol,v1,v2] = scgprops(fig);
  [x,c] = hpparam(w,beta,[],[trace,tol]);
  disp('Finished parameter problem.')
  scgset(fig, 'prevertices',x, 'const',c, 'maptype','hp2p')
  scgenable(fig,2,'on');
  
elseif strcmp(func,'d2p')		%** disk -> polygon
  [w,beta] = scgget(fig, 'vertices','angles');
  [w,beta] = scfix('d',w,beta);
  scgset(fig,'ver',w,'ang',beta)
  disp('Solving parameter problem...')
  [trace,tol,v1,v2] = scgprops(fig);
  [z,c] = dparam(w,beta,[],[trace,tol]);
  disp('Finished parameter problem.')
  scgset(fig, 'prevertices',z, 'const',c, 'maptype','d2p')
  scgenable(fig,2,'on');

elseif strcmp(func,'d2ep')		%** disk -> exterior polygon
  [w,beta] = scgget(fig, 'vertices','angles');
  [w,beta] = scfix('de',w,beta);
  scgset(fig,'ver',w,'ang',beta)
  disp('Solving parameter problem...')
  [trace,tol,v1,v2] = scgprops(fig);
  [z,c] = deparam(w,beta,[],[trace,tol]);
  disp('Finished parameter problem.')
  scgset(fig, 'prevertices',z, 'const',c, 'maptype','d2ep')
  scgenable(fig,2,'on');

elseif strcmp(func,'st2p')		%** strip -> polygon
  [w,beta] = scgget(fig, 'vertices','angles');
  disp('Use mouse to select images of left and right ends of the strip.')
  figure(gcf)
  ends = scselect(w,beta,2);
  [w,beta,ends] = scfix('st',w,beta,ends);
  scgset(fig,'ver',w,'ang',beta)
  disp('Solving parameter problem...')
  [trace,tol,v1,v2] = scgprops(fig);
  [z,c] = stparam(w,beta,ends,[],[trace,tol]);
  disp('Finished parameter problem.')
  scgset(fig, 'prevertices',z, 'const',c, 'maptype','st2p')
  scgenable(fig,2,'on');

elseif strcmp(func,'r2p')		%** rectangle -> polygon
  [w,beta] = scgget(fig, 'vertices','angles');
  disp('Use mouse to select images of rectangle corners.')
  disp('Go in counterclockwise order and select a long edge first.')
  figure(gcf)
  corner = scselect(w,beta,4);
  [w,beta,corner] = scfix('r',w,beta,corner);
  scgset(fig,'ver',w,'ang',beta)
  disp('Solving parameter problem...')
  [trace,tol,v1,v2] = scgprops(fig);
  [z,c,L] = rparam(w,beta,corner,[],[trace,tol]);
  disp('Finished parameter problem.')
  scgset(fig, 'prevertices',z, 'const',c, 'maptype','r2p')
  menus = get(fig,'userdata');
  set(menus(17),'userdata',L)
  scgenable(fig,2,'on');

elseif strcmp(func,'contin')		%** continuation
  [w,beta,z,c,maptype] = scgget(fig);
  n = length(w);
  [w,beta,idx] = modpoly(w,beta);
  if any(isnan(idx)) | any(diff([0;idx;n+1])~=1)
    fprintf('\nCannot continue after vertices have been added or deleted.\n')
    fprintf('Use direct solution instead.\n')
    scgset(fig,'ver',w,'ang',beta)
    scgenable(fig,2,'off')
    return
  end
  z0 = z;
  [trace,tol,v1,v2] = scgprops(fig);
  disp('Solving parameter problem...')
  if strcmp(maptype,'r2p')
    menus = get(fig,'userdata');
    L= get(menus(17),'userdata');
    [w,beta,z0,corners] = rcorners(w,beta,z0);
    z0 = r2strip(z0,z,L);
    [z,c,L] = rparam(w,beta,corners,z0,[trace,tol]);
    set(menus(17),'userdata',L)
  elseif strcmp(maptype,'st2p')
    ends = [find(isinf(z0)&(z0<0)),find(isinf(z0)&(z0>0))];
    [z,c] = stparam(w,beta,ends,z0,[trace,tol]);
  else
    m = size(maptypes,1);
    for j = 1:m 			% find correct prefix
      if strcmp(maptype,deblank(maptypes(j,:)))
	eval(['[z,c]=',...
		deblank(prefixes(j,:)),'param(w,beta,z0,[trace,tol]);']);
	break
      end
    end  
  end
  disp('Finished parameter problem.')
  scgset(fig,'ver',w,'ang',beta,'pre',z,'const',c,'map',maptype)
  
elseif strcmp(func,'disp') 		%** pretty print
  [w,beta,z,c,maptype] = scgget(fig);
  m = size(maptypes,1);
  if strcmp(maptype,'r2p')   % rectangle is special
    menus = get(fig,'userdata');
    L= get(menus(17),'userdata');
    rdisp(w,beta,z,c,L);
  else
    for j = 1:m 			% find correct prefix
      if strcmp(maptype,deblank(maptypes(j,:)))
	eval([deblank(prefixes(j,:)),'disp(w,beta,z,c)']);
	break
      end
    end  
  end
  
elseif strcmp(func,'plot') 		%** plot images of grid
  [w,beta,z,c,maptype] = scgget(fig,'ver','ang','pre','con','map');
  [trace,tol,v1,v2] = scgprops(fig);
  if strcmp(maptype,'r2p')
    menus = get(fig,'userdata');
    L =get(menus(17),'userdata');
    rplot(w,beta,z,c,L,v1,v2)
  else
    m = size(maptypes,1);
    for j = 1:m 			% find correct prefix
      if strcmp(maptype,deblank(maptypes(j,:)))
	eval([deblank(prefixes(j,:)),'plot(w,beta,z,c,v1,v2)']);
	break
      end
    end
  end
  disp('Finished plot.')

elseif strcmp(func,'source')		%** point source
  [w,beta,z,c,mtype] = scgget(fig);
  [trace,tol,v1,v2] = scgprops(fig);
  if strcmp(mtype,'hp2p')
    [z,c] = hp2disk(w,beta,z,c);
    mtype = 'd2p';
  end
  if strcmp(mtype,'d2p')
    ptsource(w,beta,z,c,[],v1,v2);
  else
    ptsource(w,beta,[],[],[],v1,v2);
  end
  
elseif strcmp(func,'prop')		%** Properties window
  menus = get(fig,'userdata');
  data = get(menus(4,1),'userdata');
  propfig = data(1);
  deleted = 0;
  eval('get(propfig,''pos'');','deleted=1;');
  if deleted
    screen = get(0,'screensize');
    pos = [100, screen(4)-210, 300,200];
    propfig = figure('numbertitle','off','name','SC Properties','pos',pos);
    uicontrol('style','frame','units','norm','pos',[0 0 1 1]);
    uicontrol('style','push','pos',[120,10,60,20],'string','Done',...
	'call','set(gcf,''vis'',''off'')')
    uicontrol('style','check','pos',[20,170,260,20],...
	'string','Trace parameter problem solution',...
	'call','scguicb(''pr_01'')','value',data(3));
    uicontrol('style','text','pos',[20,140,120,20],...
	'string','Error tolerace: 1e-');
    uicontrol('style','edit','pos',[143,140,20,20],...
	'string',int2str(data(4)),'call','scguicb(''pr_02'')');
    uicontrol('style','text','pos',[20,110,160,20],...
	'string','Number of curves to plot:');
    uicontrol('style','text','pos',[60,85,103,20],...
	'string','vertical/circular:');
    uicontrol('style','edit','pos',[170,85,20,20],...
	'string',int2str(data(5)),'call','scguicb(''pr_03'')');
    uicontrol('style','text','pos',[60,60,105,20],...
	'string','horizontal/radial:');
    uicontrol('style','edit','pos',[170,60,20,20],...
	'string',int2str(data(6)),'call','scguicb(''pr_04'')');
    set(menus(4,1),'userdata',[propfig,data(2:6)]);
    set(propfig,'user',[propfig,data(2:6)]);
    drawnow
  else
    set(propfig,'vis','on')
    figure(propfig)
  end
  
elseif strcmp(func(1:3),'pr_')		%** set properties
  data = get(gcf,'user');
  propfig = data(1);
  fig = data(2);
  propnum = eval(func(4:5));
  ctrl = get(propfig,'currentobject');
  if propnum==1,
    data(3) = get(ctrl,'value');
  else 
    data(2+propnum) = eval(get(ctrl,'string'));
  end
  set(propfig,'user',data);
  menus = get(fig,'user');
  set(menus(4,1),'user',data)
    
  
  
end

  


