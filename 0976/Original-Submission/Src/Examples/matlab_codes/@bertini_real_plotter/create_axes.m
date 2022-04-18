function create_axes(br_plotter)


	curr_axes = gca;

	set(curr_axes,'visible','on');
	set(curr_axes,'FontSize',br_plotter.options.fontsizes.axis);

	hold(curr_axes,'on');

	set(curr_axes,'Tag','bertini_real');

	set(curr_axes,'Position', [0.1 0.1 0.8 0.8]);
	br_plotter.axes.main = curr_axes;

end