function br_plotter = controls(br_plotter)




set_initial_visibility(br_plotter);

visibility_setup(br_plotter);

twiddle_visibility(br_plotter);


if length(br_plotter.indices)==3
    camera_setup(br_plotter);
end

end
