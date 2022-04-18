function scgset(fig,p1,v1,p2,v2,p3,v3,p4,v4,p5,v5)
%SCGSET Set data in the SCM Toolbox GUI.
%	SCGSET(FIG,'property',VALUE) sets the value of the specified
%	property to VALUE in the Schwarz-Christoffel Mapping Toolbox GUI in
%	igure FIG.  Valid properties are:
%	
%	   vertices    (polygon vertices)
%	   angles      (turning angles)
%	   prevertices (solution of parameter problem)
%	   constant    (multiplicative constant)
%	   maptype     ('hp2p', 'd2p', 'd2ep', 'st2p', 'r2p')
%	
%	Only the first three characters need be specified.  If additional
%	property name-value pairs are given, they will be set appropriately.
%	
%	SCGSET(FIG,'clear') removes all data associated with the GUI.
%
%	See also SCGGET, SCGUI.
%	
%	Written by Toby Driscoll.  Last updated 5/26/95.

maptypes = str2mat('hp2p','d2p','d2ep','st2p','r2p');

menus = get(fig,'userdata');

if (nargin==2) & strcmp(lower(p1),'clear')
  set(menus(1,1),'userdata',[]);
  return
end

if rem(nargin,2) ~= 1
  error('Wrong number of input parameters.')
end
data = get(menus(1,1),'userdata');

for k = 1:(nargin-1)/2
  prop = eval(['p',int2str(k)]);
  if ~isstr(prop)
    error('Property name expected.')
  else
    prop = lower(prop);
  end
  val = eval(['v',int2str(k)]);
  if isempty(val)
    if strcmp(prop(1:3),'con')  
      val = 0;
    elseif ~strcmp(prop(1:3),'map')
      val = zeros(size(data(:,1)));
    end
  end
      
  if strcmp(prop(1:3),'ver')
    data = [];
    data(1:length(val),1) = val(:);
  elseif strcmp(prop(1:3),'ang')
    data(1:length(val),2) = val(:);
  elseif strcmp(prop(1:3),'pre')
    data(1:length(val),3) = val(:);
  elseif strcmp(prop(1:3),'con')
    if length(val) > 1
      error('Invalid value for property ''constant''.')
    end
    data(1,4) = val;
  elseif strcmp(prop(1:3),'map')
    found = 0;
    [m,n] = size(maptypes);
    for j = 1:m
      if strcmp(lower(val),deblank(maptypes(j,:)))
	data(2,4) = j;
	found = 1;
	break
      end
    end
    if ~found & ~isempty(val)
      error(['Map type ',val,' unknown.'])
    end
  else
    error(['Property ',prop,' unknown.'])
  end
end

set(menus(1,1),'userdata',data);

