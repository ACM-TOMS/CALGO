function visibility_setup(br_plotter)





cb_params.y_start = 2;
cb_params.x = 0;
cb_params.w = 150;
cb_params.h = 15;

cb_params.curr_y = cb_params.y_start;


cb_params.y_pad = 0;
cb_params.x_pad = 0;




cb_params = make_common_checks(br_plotter,cb_params);

if br_plotter.options.render_vertices
	cb_params = make_vertex_checks(br_plotter, cb_params);
end


switch br_plotter.dimension
	
	
	case 1 % a curve
		cb_params = make_curve_checks(br_plotter, cb_params);
        

	case 2 % a surface
		if or(br_plotter.options.render_faces, br_plotter.options.render_curves)
			cb_params = make_surface_checks(br_plotter, cb_params);
		end
		%end 3d case
	otherwise
		
		
end



	
end











function cb_params = make_common_checks(br_plotter, cb_params)


h = uipanel('units','pixels','visible','on');

cb_params.curr_y = cb_params.y_start;

if ~isempty(br_plotter.legend)
	[br_plotter.checkboxes.legend, cb_params] = make_switch_checkbox('legend', 'legend', cb_params, h, br_plotter);
end

[br_plotter.checkboxes.main_axes_visibility, cb_params] = make_switch_checkbox('main axes', 'main_axes', cb_params, h, br_plotter);


if br_plotter.is_bounded == 0
	[br_plotter.checkboxes.sphere, cb_params] = make_switch_checkbox('bounding sphere', 'show_sphere', cb_params, h, br_plotter);
end


[br_plotter.checkboxes.projection_visibility, cb_params] = make_switch_checkbox('projection', 'display_projection', cb_params, h, br_plotter);



set_pos(h,cb_params);

br_plotter.panels.common_visibility = h;
cb_params.last_handle = h;
end



function cb_params = make_surface_checks(br_plotter, cb_params)
% this function build the set of checks from the bottom up.  so last checks
% are at the top of the panel.
h = uipanel('units','pixels','visible','on');
cb_params.curr_y = cb_params.y_start;



f = fieldnames(br_plotter.handles.refinements);
have_refinements = 0;
for ii = 1:length(f)
	if ~isempty(br_plotter.handles.refinements.(f{ii}))
		have_refinements = 1;
		break
	end
end










if br_plotter.options.labels
	if br_plotter.options.render_curves
		[br_plotter.checkboxes.label_critcurve, cb_params] = make_switch_checkbox('critcurve labels', 'label_critcurve', cb_params, h, br_plotter);
		[br_plotter.checkboxes.label_spherecurve, cb_params] = make_switch_checkbox('spherecurve labels', 'label_spherecurve', cb_params, h, br_plotter);




		if ~isempty(br_plotter.handles.midtext)
			new_color = get(br_plotter.handles.midtext);
			new_color = new_color.Color;
		else
			new_color = [0 0 0];
		end

		[br_plotter.checkboxes.label_midedges, cb_params] = make_switch_checkbox('mideedge labels', 'label_midedges', cb_params, h, br_plotter);
		set(br_plotter.checkboxes.label_midedges,'ForeGroundColor',new_color(1,:));



		if ~isempty(br_plotter.handles.crittext)
			new_color = get(br_plotter.handles.crittext);
			new_color = new_color.Color;
		else
			new_color = [0 0 0];
		end

		[br_plotter.checkboxes.label_critedges, cb_params] = make_switch_checkbox('criteedge labels', 'label_critedges', cb_params, h, br_plotter);
		set(br_plotter.checkboxes.label_critedges,'ForeGroundColor',new_color(1,:));



		if ~isempty(br_plotter.handles.singular_curves)
			if ~isempty(br_plotter.handles.singtext)
				new_color = get(br_plotter.handles.singtext);
				new_color = new_color.Color;
			else
				new_color = [0 0 0];
			end

			[br_plotter.checkboxes.label_singular, cb_params] = make_switch_checkbox('singedge labels', 'label_singular', cb_params, h, br_plotter);
			set(br_plotter.checkboxes.label_critedges,'ForeGroundColor',new_color(1,:));
		end

	end
	if br_plotter.options.labels
		if ~isempty(br_plotter.handles.faces)
			[br_plotter.checkboxes.label_faces, cb_params] = make_switch_checkbox('face labels', 'label_faces', cb_params, h, br_plotter);
		end
	end

end %if options.labels




cb_params.curr_y = cb_params.curr_y+10;



cb_params.x_pad = 10;

if br_plotter.options.render_curves
	[br_plotter.checkboxes.show_critcurve, cb_params] = make_switch_checkbox('crit curve', 'show_critcurve', cb_params, h, br_plotter);


	if ~br_plotter.is_bounded
		[br_plotter.checkboxes.show_spherecurve, cb_params] = make_switch_checkbox('sphere curve', 'show_spherecurve', cb_params, h, br_plotter);
	end


	[br_plotter.checkboxes.show_critslices, cb_params] = make_switch_checkbox('critslices', 'show_critslices', cb_params, h, br_plotter);

	[br_plotter.checkboxes.show_midslices, cb_params] = make_switch_checkbox('midslices', 'show_midslices', cb_params, h, br_plotter);



	if ~isempty(br_plotter.handles.singular_curves)
		[br_plotter.checkboxes.show_singular, cb_params] = make_switch_checkbox('singular curves', 'show_singular', cb_params, h, br_plotter);
	end

	cb_params.x_pad = 0;
	[br_plotter.checkboxes.surface_curve_raw, cb_params] = make_switch_checkbox('raw curves', 'raw_curves_main', cb_params, h, br_plotter);

end % re: if render_curves


if ~isempty(br_plotter.handles.faces)
	[br_plotter.checkboxes.display_faces, cb_params] = make_switch_checkbox('raw faces', 'display_faces', cb_params, h, br_plotter);
end


if or(br_plotter.options.render_faces,br_plotter.options.render_curves)
	cb_params.curr_y = cb_params.curr_y+10;
end

if have_refinements
	
	

	
	cb_params.x_pad = 10;
	
	%[checkbox_handle, cb_params] = make_switch_checkbox(switch_text, switch_name, cb_params, panel_handle, br_plotter, groupname)
	
	if br_plotter.options.render_curves
		if ~isempty(br_plotter.handles.singular_curves)
			[br_plotter.checkboxes.show_singcurve_refinement, cb_params] = make_switch_checkbox('singular curve refinements', 'show_singularcurve', cb_params, h, br_plotter,'curve_refinements');
		end

		if ~br_plotter.is_bounded
			[br_plotter.checkboxes.show_spherecurve_refinement, cb_params] = make_switch_checkbox('sphere curve refinements', 'show_spherecurve', cb_params, h, br_plotter,'curve_refinements');
		end
		[br_plotter.checkboxes.show_critcurve_refinement, cb_params] = make_switch_checkbox('crit curve refinements', 'show_critcurve', cb_params, h, br_plotter,'curve_refinements');
		[br_plotter.checkboxes.show_critslice_refinement, cb_params] = make_switch_checkbox('crit slice refinements', 'show_critslices', cb_params, h, br_plotter,'curve_refinements');
		[br_plotter.checkboxes.show_midslice_refinement, cb_params] = make_switch_checkbox('mid slice refinements', 'show_midslices', cb_params, h, br_plotter,'curve_refinements');

		cb_params.x_pad = 0;

		[br_plotter.checkboxes.curve_refinements, cb_params] = make_switch_checkbox('curve refinements', 'main', cb_params, h, br_plotter,'curve_refinements');
	end%re: if render_curves
	
	if ~isempty(br_plotter.handles.faces)
		[br_plotter.checkboxes.display_face_samples, cb_params] = make_switch_checkbox('face samples', 'display_face_samples', cb_params, h, br_plotter);
	end
end




set_pos(h,cb_params);

br_plotter.panels.surface = h;
cb_params.last_handle = h;
end





function cb_params = make_vertex_checks(br_plotter, cb_params)

h = uipanel('units','pixels','visible','on');

cb_params.curr_y = cb_params.y_start;



[br_plotter.checkboxes.vertex_marker_visibility, cb_params] = make_switch_checkbox('vertex markers', 'display_vertices', cb_params, h, br_plotter);

if br_plotter.options.labels
	[br_plotter.checkboxes.vertex_label_visibility, cb_params] = make_switch_checkbox('vertex labels', 'label_vertices', cb_params, h, br_plotter);
end %re if labels
	f = br_plotter.legend.vertices.types;
	for ii = 1:length(f)

		[br_plotter.checkboxes.vertex_set_visibility(ii), cb_params] = make_switch_checkbox(f{ii}, f{ii}, cb_params, h, br_plotter,'vertex_set');
		
		
		new_color = get(br_plotter.handles.vertices.(f{ii}));
		new_color = new_color.Color;

		set(br_plotter.checkboxes.vertex_set_visibility(ii),'ForeGroundColor',new_color(1,:))



	end


set_pos(h,cb_params);

br_plotter.panels.vertex = h;
cb_params.last_handle = h;
end









function cb_params = make_curve_checks(br_plotter, cb_params)

h = uipanel('units','pixels','visible','on');





[br_plotter.checkboxes.critcurve_visibility, cb_params] = make_switch_checkbox('edges', 'show_edges', cb_params, h, br_plotter);
if ~isempty(br_plotter.handles.sample_edges)
	[br_plotter.checkboxes.critcurve_samples, cb_params] = make_switch_checkbox('refinement', 'show_curve_samples', cb_params, h, br_plotter);
end




set_pos(h,cb_params);

cb_params.last_handle = h;
br_plotter.panels.curve = h;

end










function set_pos(h, cb_params)


if isfield(cb_params,'last_handle')
	prev_pos = get(cb_params.last_handle,'Position');
	start_y = prev_pos(2) + prev_pos(4) + 5;
else
	start_y = 5;
end
	
pos = [5 start_y cb_params.w+3*cb_params.x_pad cb_params.curr_y+cb_params.y_pad+cb_params.y_start+5];


set(h,'Position',pos);


end



function [checkbox_handle, cb_params] = make_switch_checkbox(switch_text, switch_name, cb_params, panel_handle, br_plotter, groupname)


position =[cb_params.x_pad cb_params.curr_y cb_params.w-5-3*cb_params.x_pad cb_params.h];

if nargin == 6
	
checkbox_handle = uicontrol(...
				'style','checkbox',...
				'units','pixels',...
				'position',position,...
				'String',switch_text,...
				'Value',br_plotter.switches.(groupname).(switch_name),...
				'callback',@(src,event)flip_switch(br_plotter,src,event,switch_name,groupname),...
				'Parent',panel_handle);
else
	
checkbox_handle = uicontrol(...
				'style','checkbox',...
				'units','pixels',...
				'position',position,...
				'String',switch_text,...
				'Value',br_plotter.switches.(switch_name),...
				'callback',@(src,event)flip_switch(br_plotter,src,event,switch_name),...
				'Parent',panel_handle);
end


			
cb_params.curr_y = cb_params.curr_y+cb_params.h+cb_params.y_pad;


end




