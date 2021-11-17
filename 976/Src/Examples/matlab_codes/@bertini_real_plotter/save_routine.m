function br_plotter = save_routine(br_plotter,varargin) % additional arguments go afer the first two, which must be put here ,even if they are not used.
%perhaps get info from other calls here?

f = fieldnames(br_plotter.panels);

for ii = 1:length(f)
	set( findall(br_plotter.panels.(f{ii}), '-property', 'visible'), 'visible', 'off')
	set(br_plotter.panels.(f{ii}),'visible','off');
	
end

try
	if and(br_plotter.options.render_faces,or(br_plotter.switches.display_faces == 1,br_plotter.switches.display_face_samples == 1))
		br_plotter.options.format = 'png';
		br_plotter.options.format_flag = 'png';
		render_into_file(br_plotter.options,'-r300');
	else
		br_plotter.options.format = 'eps';
		br_plotter.options.format_flag = 'psc2';
		render_into_file(br_plotter.options);
	end
catch
	warning('unable to complete saving for some reason');
end

for ii = 1:length(f)
	set(br_plotter.panels.(f{ii}),'visible','on');
	set( findall(br_plotter.panels.(f{ii}), '-property', 'visible'), 'visible', 'on')
end

end
