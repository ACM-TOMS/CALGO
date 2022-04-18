function center_camera_on_selected_point(br_plotter,source, event)


curr_axes = br_plotter.axes.main;

dcm_obj = datacursormode(br_plotter.figures.main);
set(dcm_obj,'DisplayStyle','datatip',...
	'SnapToDataVertex','on','Enable','on');

set(curr_axes,'CameraPositionMode','manual');

curr_pos = get(curr_axes,'CameraPosition');


cursor_info = getCursorInfo(dcm_obj);


if isfield(cursor_info,'Position')
	

	br_plotter.scene.target = cursor_info.Position;
	set(curr_axes,'CameraTarget',cursor_info.Position);
	set(curr_axes,'CameraPosition',curr_pos);
else
	helpdlg('Select a datapoint on which to center, then try again',...
        'Center Point Selection');
end






	
end
