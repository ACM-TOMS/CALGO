function br_plotter = camera_setup(br_plotter)

cameratoolbar

curr_axes = br_plotter.axes.main;


set(curr_axes,'CameraTargetMode','manual');
br_plotter.scene.target = real(br_plotter.BRinfo.center(br_plotter.indices));
set(curr_axes,'CameraTarget',br_plotter.scene.target);

new_campos = [get(curr_axes,'XLim');get(curr_axes,'YLim');get(curr_axes,'ZLim')];

br_plotter.scene.campos = 5*new_campos([2 3 6]);

set(curr_axes,'CameraPosition',br_plotter.scene.campos);
set(curr_axes,'CameraViewAngle',10);

rotate3d off

end
