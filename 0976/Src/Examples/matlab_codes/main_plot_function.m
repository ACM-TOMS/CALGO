function h = main_plot_function(data,ind, curr_axis)
%
switch length(ind)
	case {2}
		h = plot(data(:,ind(1)),data(:,ind(2)),'-','Parent',curr_axis);
	case {3}
		h = plot3(data(:,ind(1)),data(:,ind(2)),data(:,ind(3)),'-','Parent',curr_axis);
	otherwise
			
end
	
end


