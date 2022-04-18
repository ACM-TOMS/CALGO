function br_plotter = change_alpha(br_plotter,source,event)
			
			
			prompt = {'New face alpha:','New refinement alpha:'};
			dlg_title = 'Alpha';
			num_lines = 2;
			def = {num2str(br_plotter.options.face_alpha),...
				   num2str(br_plotter.options.sample_alpha)};

			output = inputdlg(prompt,dlg_title,num_lines,def);

			if isempty(output)
				return;
			end
			
			new_face_alpha = str2num(output{1});
			new_sample_alpha = str2num(output{2});
			
			if br_plotter.options.face_alpha ~= new_face_alpha
				set(br_plotter.handles.faces(:),'FaceAlpha',new_face_alpha);
				br_plotter.options.face_alpha = new_face_alpha;
			end
			
			if br_plotter.options.sample_alpha ~= new_sample_alpha
				set(br_plotter.handles.surface_samples(:),'FaceAlpha',new_sample_alpha);
				br_plotter.options.sample_alpha = new_sample_alpha;
			end
			
			set(br_plotter.buttons.alpha,'String',sprintf('Alpha (%i/%i)',mod(100*br_plotter.options.face_alpha,100),mod(100*br_plotter.options.sample_alpha,100)));
			

end