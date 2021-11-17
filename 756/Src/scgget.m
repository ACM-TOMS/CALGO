function [val1,val2,val3,val4,val5] = scgget(fig,p1,p2,p3,p4,p5)
%SCGGET Get data from the SC Toolbox GUI.
%	SCGGET(FIG,'property') returns the value of the specified property
%	associated with the Schwarz-Christoffel Toolbox GUI in
%	figure FIG.  The properties are:
%	
%	   vertices    (polygon vertices)
%	   angles      (turning angles)
%	   prevertices (solution of parameter problem)
%	   constant    (multiplicative constant)
%	   maptype     ('hp2p', 'd2p', 'd2ep', 'st2p', 'r2p')
%	
%	Only the first three characters need be specified.  If additional
%	property names are given, values will be returned in additional
%	output arguments in the same order.
%	
%	SCGGET(FIG) is shorthand for SCGGET(FIG,'ver','ang','pre','con',...
%	'map').
%
%       See also SCGSET, SCGUI.
%
%	Written by Toby Driscoll.  Last updated 5/24/95.


if nargin==1
  p1 = 'ver';
  p2 = 'ang';
  p3 = 'pre';
  p4 = 'con';
  p5 = 'map';
  nargin = 6;
elseif nargout~=(nargin-1) & nargout~=0
  error('Incorrect number of output parameters.')
end

maptypes = str2mat('hp2p','d2p','d2ep','st2p','r2p');

menus = get(fig,'userdata');
data = get(menus(1,1),'userdata');
[n,p] = size(data);
if p==0
  return
end

w = data(:,1);
if ~any(w), w = []; end;
beta = [];
z = [];
c = [];
mapnum = 0;
if p > 1
  beta = data(:,2);
  if ~any(beta), beta = []; end
  if p > 2
    z = data(:,3);
    if ~any(z), z = []; end
    if p > 3
      c = data(1,4);
      if ~any(c), c = []; end
      mapnum = data(2,4);
    end
  end
end

for k = 1:(nargin-1)
  prop = eval(['p',int2str(k)]);
  if strcmp(lower(prop(1:3)),'ver')
    val = w;
  elseif strcmp(lower(prop(1:3)),'ang')
    val = beta;
  elseif strcmp(lower(prop(1:3)),'pre')
    val = z;
  elseif strcmp(lower(prop(1:3)),'con')
    val = c;
  elseif strcmp(lower(prop(1:3)),'map')
    if ~mapnum
      val = [];
    else
      val = deblank(maptypes(mapnum,:));
    end
  end
  eval(['val',int2str(k),' = val;']);
end



