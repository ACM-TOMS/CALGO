% the callback function for resizing the ui.
%
% keeps things looking nice as the window is resized.
%
% dani brake
% danielthebrake@gmail.com
% 2014-2016

function resizeui(br_plotter,srcHandle,eventData,varargin)
	br_plotter.window = get(br_plotter.figures.main,'Position');
	w = br_plotter.window;
	if ~isempty(br_plotter.panels)

		p = get(br_plotter.panels.buttons,'Position');
		set(br_plotter.panels.buttons,'position',[5    w(4)-p(4)-5     p(3)    p(4)]);


		c = get(br_plotter.panels.common_visibility,'Position');
		total_vertical_size = c(4);
		if isfield(br_plotter.panels,'vertex')
			v = get(br_plotter.panels.vertex,'Position');
			total_vertical_size = total_vertical_size+v(4);
		end

		if isfield(br_plotter.panels,'surface')
			s = get(br_plotter.panels.surface,'Position');
			total_vertical_size = total_vertical_size+s(4);
		end

		if isfield(br_plotter.panels,'curve')
			cu = get(br_plotter.panels.curve,'Position');
			total_vertical_size = total_vertical_size+cu(4);
		end


		if total_vertical_size+30 > w(4)-p(4)
			P = w(3)-c(3)-5;
		else
			P = 5;
		end

		set(br_plotter.panels.common_visibility,'position',[P 5 c(3)    c(4)]);
		Q = 10+c(4);
		if isfield(br_plotter.panels,'vertex')
			set(br_plotter.panels.vertex,'position',[P Q v(3)    v(4)]);
			Q = Q+5+v(4);
		end

		if isfield(br_plotter.panels,'surface')
			set(br_plotter.panels.surface,'position',[P Q s(3)    s(4)]);
			Q = Q+5+s(4);
		end

		if isfield(br_plotter.panels,'curve')
			set(br_plotter.panels.curve,'position',[P Q cu(3)    cu(4)]);
			Q = Q+5+cu(4);
		end				

	end
end	