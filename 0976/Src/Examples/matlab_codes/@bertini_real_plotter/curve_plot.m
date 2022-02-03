%curve specific code

function br_plotter = curve_plot(br_plotter)
%

create_axes(br_plotter);

label_axes(br_plotter);

% actually plot
plot_projection(br_plotter);


plot_vertices(br_plotter);


adjust_axes(br_plotter);


br_plotter.handles.sample_edges = [];
br_plotter.handles.edges = [];

plot_edge_points(br_plotter);


if ~isempty(br_plotter.BRinfo.sampler_data)
    
    if isfield(br_plotter.options,'use_fixed_linestyle')
        style = br_plotter.options.linestyle;
    else
        style = '-';
    end


	br_plotter.handles.sample_edges = plot_curve_samples(br_plotter,br_plotter.BRinfo.sampler_data,style);
end

sphere_plot(br_plotter);

end


























