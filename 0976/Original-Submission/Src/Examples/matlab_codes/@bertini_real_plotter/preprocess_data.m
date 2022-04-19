function preprocess_data(br_plotter)

	for ii = 1:br_plotter.BRinfo.num_vertices
		temp = br_plotter.options.custom_projection(br_plotter.BRinfo.vertices(ii).point);
		br_plotter.BRinfo.vertices(ii).point = temp;
	end

	br_plotter.BRinfo.var_names = {};
	for ii = 1:length(br_plotter.BRinfo.vertices(1).point)
		br_plotter.BRinfo.var_names{ii} = ['proj_' num2str(ii)];
	end

	
	br_plotter.BRinfo.num_variables = length(br_plotter.BRinfo.vertices(1).point)+1;

end