%  surface specific code



function br_plotter = surface_plot(br_plotter)

create_axes(br_plotter);

label_axes(br_plotter);

sphere_plot(br_plotter);

plot_vertices(br_plotter);

adjust_axes(br_plotter);


if br_plotter.options.render_curves
	plot_surface_edges(br_plotter);
end

if br_plotter.options.render_faces
	plot_faces(br_plotter);
	plot_surface_samples(br_plotter);
end

plot_projection(br_plotter);




sync_axes(br_plotter);


end












