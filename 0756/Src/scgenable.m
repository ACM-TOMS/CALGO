function scgenable(fig,type,action)
%SCGENABLE Enables/disables menus in the SCM Toolbox GUI.
%	Menu items created by SCGUI are at times disabled, when they
%	appear to have no current meaning.  However, if you use SCGSET,
%	conditions may change.  SCGENABLE(FIG,TYPE,ACTION) will change
%	the status of a class(es) of menus.  TYPE is an integer.
%	Generally, TYPE=1 is on when a polygon is known to the GUI, and
%	TYPE=2 is on when a parameter problem solution is known.  ACTION
%	is either 'on' or 'off'.
%
%	Written by Toby Driscoll.  Last updated 5/23/95.

%       This whole mechanism could be a lot friendlier.

menus = get(fig,'userdata');
for i = 1:length(type);
  for j = find(menus(:,2)==type(i))
    set(menus(j,1), 'enable',action)
  end
end
