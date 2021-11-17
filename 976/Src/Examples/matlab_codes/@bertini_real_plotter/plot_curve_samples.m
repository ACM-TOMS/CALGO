function handles = plot_curve_samples(br_plotter,sampler_data,style,desiredcolor)

if nargin == 4
	repeat_colors = true;
else
	repeat_colors = false;
end

ind = br_plotter.indices;



num_non_degen = 0;

indicator = zeros(size(sampler_data.sample_sizes,1),1);
for ii = 1:size(sampler_data.sample_sizes,1)
	
	left_point = br_plotter.BRinfo.vertices(sampler_data.edge(ii).samples(1)+1).point(1:br_plotter.BRinfo.num_variables-1);
	right_point = br_plotter.BRinfo.vertices(sampler_data.edge(ii).samples(end)+1).point(1:br_plotter.BRinfo.num_variables-1);
	
	if norm(left_point-right_point)<1e-8
		continue;
	else
		indicator(ii) = 1;
		num_non_degen = num_non_degen+1;
	end
	
end

if repeat_colors
	colors = repmat(desiredcolor,[num_non_degen 1]);
else
	colors = br_plotter.options.colormap(num_non_degen);
end



curr_axis = br_plotter.axes.main;
nondegen_edge_ind = 1;

handles = zeros(num_non_degen,1);

for ii = 1:length(sampler_data.edge)
	
	curr_samples = sampler_data.edge(ii).samples;
	
	if indicator(ii)==0
       continue; 
	end
    
	plotme = br_plotter.fv.vertices(curr_samples+1,:);
	

	h = main_plot_function(plotme,1:length(ind),curr_axis);

	
	set(h,'Color',colors(nondegen_edge_ind,:));
	set(h,'LineWidth',br_plotter.options.line_thickness);
    set(h,'LineStyle',style);
	
	
    handles(nondegen_edge_ind) = h;
	
	
	nondegen_edge_ind = nondegen_edge_ind+1;
	
end




end













