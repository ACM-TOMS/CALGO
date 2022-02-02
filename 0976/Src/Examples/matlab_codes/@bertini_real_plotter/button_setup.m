function button_setup(br_plotter)



h = uipanel('units','pixels','visible','on');
br_plotter.panels.buttons = h;


height = 20;
width = 100;

vert_spacing = 5;
horiz_pad = 5;
curr_y = 0;



[curr_y, br_plotter.buttons.save] = make_button('Save', @br_plotter.save_routine, horiz_pad, vert_spacing, curr_y, br_plotter.panels.buttons, height, width);

[curr_y, br_plotter.buttons.load] = make_button('Load & Render', {@br_plotter.load_and_render}, horiz_pad, vert_spacing, curr_y, br_plotter.panels.buttons, height, width);


[curr_y, br_plotter.buttons.center] = make_button('Center Camera', {@br_plotter.center_camera_on_selected_point}, horiz_pad, vert_spacing, curr_y, br_plotter.panels.buttons, height, width);


[curr_y, br_plotter.buttons.leap] = make_button('Leap Motion', {@leap_figure_axes_control}, horiz_pad, vert_spacing, curr_y, br_plotter.panels.buttons, height, width);


if br_plotter.options.render_faces
	[curr_y, br_plotter.buttons.stl] = make_button('Save to STL', {@fv_to_stl,br_plotter.fv}, horiz_pad, vert_spacing, curr_y, br_plotter.panels.buttons, height, width);
end

[curr_y, br_plotter.buttons.fontsize] = make_button('FontSize', {@br_plotter.change_text_size}, horiz_pad, vert_spacing, curr_y, br_plotter.panels.buttons, height, width);


buttontext = sprintf('MarkerSize %i',br_plotter.options.markersize);
[curr_y, br_plotter.buttons.markersize] = make_button(buttontext, {@br_plotter.change_markersize}, horiz_pad, vert_spacing, curr_y, br_plotter.panels.buttons, height, width);




if br_plotter.dimension == 2
	buttontext = sprintf('Alpha (%i/%i)',mod(100*br_plotter.options.face_alpha,100),mod(100*br_plotter.options.sample_alpha,100));
	[curr_y, br_plotter.buttons.alpha] = make_button(buttontext, {@br_plotter.change_alpha}, horiz_pad, vert_spacing, curr_y, br_plotter.panels.buttons, height, width);
end

vertical_size = br_plotter.window(4);
set(br_plotter.panels.buttons,'position',[horiz_pad    vertical_size-curr_y-vert_spacing-2     width+2*horiz_pad+2    curr_y+vert_spacing+2]);


end





function [new_y, button_handle] = make_button(text, callback_function, horiz_pad, vert_spacing, curr_y, panel_handle, height, width)

button_handle = uicontrol('Style', 'pushbutton', 'String', text,...
	'Position', [horiz_pad curr_y+vert_spacing width height],...
	'Callback', callback_function, 'Parent',panel_handle); 
new_y = curr_y+height+vert_spacing;


end
