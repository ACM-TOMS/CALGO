function br_plotter = twiddle_visibility(br_plotter)


common_twiddle(br_plotter);

vertex_twiddle(br_plotter);





switch br_plotter.dimension
	case 1
		
		
		curve_twiddle(br_plotter);
		
		
	case 2
		
		
		surface_twiddle(br_plotter);
		
	otherwise
		
end



end




function [] = surface_twiddle(br_plotter)

surface_subcurve_twiddle(br_plotter);


surface_face_twiddle(br_plotter);


surface_refinement_twiddle(br_plotter);






end





function [] = surface_refinement_twiddle(br_plotter)


%faces
if br_plotter.switches.display_face_samples == 0
	for ii = 1:length(br_plotter.handles.surface_samples)
		set(br_plotter.handles.surface_samples(ii),'visible','off');
	end
else
	for ii = 1:length(br_plotter.handles.surface_samples)
		set(br_plotter.handles.surface_samples(ii),'visible','on');
	end
end



if br_plotter.switches.curve_refinements.main
	
	if br_plotter.switches.curve_refinements.show_critslices
		set(br_plotter.handles.refinements.critslice(:),'visible','on');
	else
		set(br_plotter.handles.refinements.critslice(:),'visible','off');
	end
	
	
	if br_plotter.switches.curve_refinements.show_midslices
		set(br_plotter.handles.refinements.midslice(:),'visible','on');
	else
		set(br_plotter.handles.refinements.midslice(:),'visible','off');
	end
	
	
	if br_plotter.switches.curve_refinements.show_spherecurve
		set(br_plotter.handles.refinements.spherecurve,'visible','on');
	else
		set(br_plotter.handles.refinements.spherecurve,'visible','off');
	end
	
	
	if br_plotter.switches.curve_refinements.show_critcurve
		set(br_plotter.handles.refinements.critcurve,'visible','on');
	else
		set(br_plotter.handles.refinements.critcurve,'visible','off');
	end
	
	
	
	if br_plotter.switches.curve_refinements.show_singularcurve
		set(br_plotter.handles.refinements.singularcurve(:),'visible','on');
	else
		set(br_plotter.handles.refinements.singularcurve(:),'visible','off');
	end
	
	
	
else
	
	
	set(br_plotter.handles.refinements.critslice(:),'visible','off');
	set(br_plotter.handles.refinements.midslice(:),'visible','off');
	set(br_plotter.handles.refinements.spherecurve,'visible','off');
	set(br_plotter.handles.refinements.critcurve,'visible','off');
	set(br_plotter.handles.refinements.singularcurve(:),'visible','off');
	
	
	
end


end





function [] = surface_face_twiddle(br_plotter)


if br_plotter.options.labels
	if br_plotter.switches.label_faces == 0
		set(br_plotter.handles.face_labels,'visible','off');
	else
		set(br_plotter.handles.face_labels,'visible','on');
	end
end



%faces
if br_plotter.switches.display_faces == 0
	for ii = 1:length(br_plotter.handles.faces)
		set(br_plotter.handles.faces(ii),'visible','off');
	end
else
	for ii = 1:length(br_plotter.handles.faces)
		set(br_plotter.handles.faces(ii),'visible','on');
	end
end


end


function [] = surface_subcurve_twiddle(br_plotter)

if br_plotter.switches.raw_curves_main
	if br_plotter.switches.show_critcurve == 0
		set(br_plotter.handles.critcurve(:),'visible','off')
		if br_plotter.options.labels
			set(br_plotter.handles.critcurve_labels,'visible','off');
		end
	else
		set(br_plotter.handles.critcurve(:),'visible','on')
		if br_plotter.options.labels
			if (br_plotter.switches.label_critcurve == 0)
				set(br_plotter.handles.critcurve_labels,'visible','off');
			else
				set(br_plotter.handles.critcurve_labels,'visible','on');
			end
		end
	end



	if br_plotter.switches.show_spherecurve == 0
		set(br_plotter.handles.spherecurve(:),'visible','off')
		if br_plotter.options.labels
			set(br_plotter.handles.spherecurve_labels,'visible','off');
		end
	else
		set(br_plotter.handles.spherecurve(:),'visible','on')
		if br_plotter.options.labels
			if (br_plotter.switches.label_spherecurve == 0)
				set(br_plotter.handles.spherecurve_labels,'visible','off');
			else
				set(br_plotter.handles.spherecurve_labels,'visible','on');
			end
		end
	end

	if br_plotter.switches.show_critslices == 0
		set(br_plotter.handles.critslices(:),'visible','off')
		if br_plotter.options.labels
			set(br_plotter.handles.crittext,'visible','off');
		end
	else

		set(br_plotter.handles.critslices(:),'visible','on')
		if br_plotter.options.labels
			if (br_plotter.switches.label_critedges == 0)
				set(br_plotter.handles.crittext,'visible','off');
			else
				set(br_plotter.handles.crittext,'visible','on');
			end
		end
	end



	if br_plotter.switches.show_singular == 0
		set(br_plotter.handles.singular_curves(:),'visible','off')
		if br_plotter.options.labels
			set(br_plotter.handles.singtext,'visible','off');
		end
	else
		set(br_plotter.handles.singular_curves(:),'visible','on')
		if br_plotter.options.labels
			if (br_plotter.switches.label_singular == 0)
				set(br_plotter.handles.singtext,'visible','off');
			else
				set(br_plotter.handles.singtext,'visible','on');
			end
		end


	end


	if br_plotter.switches.show_midslices == 0
		set(br_plotter.handles.midslices(:),'visible','off')
		if br_plotter.options.labels
			set(br_plotter.handles.midtext,'visible','off');
		end
	else
		set(br_plotter.handles.midslices(:),'visible','on')
		if br_plotter.options.labels
			if (br_plotter.switches.label_midedges == 0)
				set(br_plotter.handles.midtext,'visible','off');
			else
				set(br_plotter.handles.midtext,'visible','on');
			end
		end


	end
	
else % turn off all curves and labels
	set(br_plotter.handles.critcurve(:),'visible','off');
	set(br_plotter.handles.critcurve_labels,'visible','off');
	
	set(br_plotter.handles.spherecurve(:),'visible','off');
	set(br_plotter.handles.spherecurve_labels,'visible','off');
	
	set(br_plotter.handles.critslices(:),'visible','off')
	set(br_plotter.handles.crittext,'visible','off');
	
	set(br_plotter.handles.midslices(:),'visible','off')
	set(br_plotter.handles.midtext,'visible','off');
	
	set(br_plotter.handles.singular_curves(:),'visible','off')
	set(br_plotter.handles.singtext,'visible','off');
end  % re main switch for raw curves



end


function [] = curve_twiddle(br_plotter)

if br_plotter.switches.show_edges == 0
	set(br_plotter.handles.edges(:),'visible','off')
else
	set(br_plotter.handles.edges(:),'visible','on')
end

if br_plotter.switches.show_curve_samples == 0
	set(br_plotter.handles.sample_edges(:),'visible','off')
else
	set(br_plotter.handles.sample_edges(:),'visible','on')
end

end





function [] = common_twiddle(br_plotter)

if ~isempty(br_plotter.legend)
	if br_plotter.switches.legend == 0
		set(br_plotter.handles.legend,'visible','off');
	else
		set(br_plotter.handles.legend,'visible','on');
	end
end

if br_plotter.switches.main_axes == 0
	set(br_plotter.axes.main,'visible','off');
else
	set(br_plotter.axes.main,'visible','on');
end


if br_plotter.switches.display_projection ==0
	set(br_plotter.handles.projection(:),'visible','off');
else
	set(br_plotter.handles.projection(:),'visible','on');
end






if br_plotter.switches.show_sphere
	set(br_plotter.handles.sphere,'visible','on');
else
	set(br_plotter.handles.sphere,'visible','off');
end

end



function [] = vertex_twiddle(br_plotter)
%vertices
if br_plotter.options.render_vertices
	
	if br_plotter.switches.display_vertices == 0
		f = fieldnames(br_plotter.switches.vertex_set);
		for ii = 1:length(f)
			set(br_plotter.handles.vertices.(f{ii}),'visible','off');
			if br_plotter.options.labels
				set(br_plotter.handles.vertex_text.(f{ii}),'visible','off');
			end
		end
		
	else
		
		f = fieldnames(br_plotter.switches.vertex_set);
		for ii = 1:length(f)
			if br_plotter.switches.vertex_set.(f{ii}) == 0
				set(br_plotter.handles.vertices.(f{ii}),'visible','off')
			else
				set(br_plotter.handles.vertices.(f{ii}),'visible','on')
			end
		end
		if br_plotter.options.labels
			for ii = 1:length(f)
				
				if br_plotter.switches.label_vertices == 0;
					
					set(br_plotter.handles.vertex_text.(f{ii}),'visible','off')
					
				else
					if br_plotter.switches.vertex_set.(f{ii}) == 0
						set(br_plotter.handles.vertex_text.(f{ii}),'visible','off');
					else
						set(br_plotter.handles.vertex_text.(f{ii}),'visible','on');
					end
				end
			end
		end
	end
	
end

end
