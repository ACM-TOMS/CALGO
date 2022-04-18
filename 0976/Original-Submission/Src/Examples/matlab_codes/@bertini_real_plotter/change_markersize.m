function change_markersize(br_plotter,src,hndl)



prompt = {'Enter new marker size:'};
			dlg_title = 'MarkerSize';
			num_lines = 1;
			default = {num2str(br_plotter.options.markersize)};

			output = inputdlg(prompt,dlg_title,num_lines,default);

			if isempty(output)
				return;
			end
			
			new_size = str2num(output{1});
			
			if new_size~=br_plotter.options.markersize
				
				f = fieldnames(br_plotter.handles.vertices);
				for ii = 1:length(f)
					set(br_plotter.handles.vertices.(f{ii}),'MarkerSize',new_size);
				end
				
				set(br_plotter.buttons.markersize,'String',sprintf('MarkerSize %i',new_size));
				br_plotter.options.markersize = new_size;
			end
			
			
			
end