function flip_switch(br_plotter,srcHandle,eventData,varargin)

switch_name = varargin{1};

if length(varargin)==2
	category_name = varargin{2};
	br_plotter.switches.(category_name).(switch_name) = mod(br_plotter.switches.(category_name).(switch_name)+1,2);
else
	br_plotter.switches.(switch_name) = mod(br_plotter.switches.(switch_name)+1,2);
end
	
	




% if length(varargin)==3
% 	
% elseif length(varargin)>3
% 	
% end




twiddle_visibility(br_plotter);
end
