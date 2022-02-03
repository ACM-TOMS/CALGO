%
% a class for plotting the data from a bertini_real run of any computable
% dimension
%
% example invokation: 
% bertini_real_plotter('autosave',false,'vertices',false,'linestyle','-','colormap',@jet,'colorfn',@norm,'num_colors',200,'curves',false)
%
%
% options: 
%	'autosave'          - bool [false]
%	'vertices', 'vert'  - bool [false]  for large samplings, this causes
%							rendering to be very slow.
%	'filename', 'file'   - string [BRinfo*.mat]
%	'proj'              - handle to function.  no default
%	'mono', 'monocolor' - color or RGB triple.  not on by default
%	'labels'            - bool [true]
%  	'colormap'          - handle to colormap generating function [@jet]
%  	'linestyle'         - string, sets all curves to have this style. 
%                            no default, curves have style by type.
%	'curves', 'curve'   - bool [true]
%	'faces'             - bool [true]
%
%   'colorfn'           - handle to function of x, for generating color
%                            data.  no default value.  if this is not
%                            specified, then the colors for the faces
%                            correspond to the entire face, in the order
%                            computed.  An example of using this colorfn
%                            would be to pass a handle to a function
%                            computing the distance between x and the
%                            origin, perhaps.
%
%    'num_colors'        - integer [64] the number of colors used in the
%							colormap, particularly used when you use the 
%							'colorfn' option, to
%							specify a function used for coloring the surface.
%							also, 'numcolors'
%
%
%
%
%
% dani brake
% danielthebrake@gmail.com
% university of notre dame
% applied and computational mathematics and statistics
% 2014, 2015, 2016



classdef bertini_real_plotter < handle
	
	
	properties
		BRinfo = [];
% 		scene = scene_manipulator();
		
		window = [20 20 1024 768];
		figures = []
		axes = [];
		handles = [];
		legend = [];
		filename = [];
		
		panels = [];
		buttons = [];
		switches = [];
		checkboxes = [];
		
		scene = [];
		
		dimension = -1;
		
		indices = [];
		
		options = [];
		
		is_bounded = [];
		fv = [];
	end
	
	methods
		
		
		function br_plotter = bertini_real_plotter(varargin)
			initialize(br_plotter);
			set_default_options(br_plotter);
			set_options(br_plotter,varargin);
			set_filename(br_plotter);
			load(br_plotter);
			plot(br_plotter);
		end %re: bertini_real_plotter() constructor
		
		function initialize(br_plotter)
			initialize_handles_surface(br_plotter);
		end
		
		initialize_handles_surface(br_plotter);
			
			
		function load(br_plotter)
			load_data(br_plotter);
			
			if br_plotter.BRinfo.num_vertices==0
				warning('your decomposition contains 0 vertices.  The real part appears to be empty.  Plotting will now terminate.')
				display(br_plotter.BRinfo);
				return;
			end
			
			if ~isfield(br_plotter.BRinfo,'run_metadata')
				br_plotter.BRinfo.run_metadata.version.number = 103;
			end
			
			if br_plotter.options.use_custom_projection
				preprocess_data(br_plotter);
			end
			
			get_indices(br_plotter);
		end
		
		function set_default_options(br_plotter)
			br_plotter.options.use_custom_projection = false;
			br_plotter.options.markersize = 10;
			br_plotter.options.sample_alpha = 1;
			br_plotter.options.face_alpha = 1;
			br_plotter.options.edge_alpha = 0.4;
			br_plotter.options.fontsizes.legend = 12;
			br_plotter.options.fontsizes.labels = 16;
			br_plotter.options.fontsizes.axis = 20;
			br_plotter.options.line_thickness = 6;
			br_plotter.options.autosave = false;
			
			br_plotter.options.labels = true;
			br_plotter.options.monocolor = false;
			
			br_plotter.options.render_vertices = false;
			br_plotter.options.render_curves = true;
			br_plotter.options.render_faces = true;
			
			br_plotter.options.use_colorfn = false;
			
			br_plotter.options.num_colors = 64;
            if isempty(which('parula'))
                br_plotter.options.colormap = @jet;
            else
                br_plotter.options.colormap = @parula;
            end
		end
		
		%parses the command line options fed into the constructor.
		function set_options(br_plotter,command_line_options)
			
			
			if mod(length(command_line_options),2)~=0
				error('must have option-value pairs');
			end
			
			
			for ii = 1:2:length(command_line_options)-1
				
				switch lower(command_line_options{ii})
					case 'autosave'
						
						tentative_arg = command_line_options{ii+1};
						
						if ischar(tentative_arg)
							switch tentative_arg
								case {'y','yes','true'}
									br_plotter.options.autosave = true;
								case {'n','no','false'}
									br_plotter.options.autosave = false;
								otherwise
									error('bad option %s for autosave',tentative_arg);
							end
							
						else
							if tentative_arg==1
								br_plotter.options.autosave = true;
							elseif tentative_arg==0
								br_plotter.options.autosave = false;
							else
								error('bad option %f for autosave',tentative_arg);
							end
						end
					
					case {'curves','curve'}
						tentative_arg = command_line_options{ii+1};
						
						if ischar(tentative_arg)
							switch tentative_arg
								case {'y','yes','true'}
									br_plotter.options.render_curves = true;
								case {'n','no','none','false'}
									br_plotter.options.render_curves = false;
								otherwise
									error('bad option %s for curves',tentative_arg);
							end
							
						else
							if tentative_arg==1
								br_plotter.options.render_curves = true;
							elseif tentative_arg==0
								br_plotter.options.render_curves = false;
							else
								error('bad option %f for curves',tentative_arg);
							end
						end
						
						
					case {'vertices','vert'}
						
						tentative_arg = command_line_options{ii+1};
						
						if ischar(tentative_arg)
							switch tentative_arg
								case {'y','yes','true'}
									br_plotter.options.render_vertices = true;
								case {'n','no','none','false'}
									br_plotter.options.render_vertices = false;
								otherwise
									error('bad option %s for vertices',tentative_arg);
							end
							
						else
							if tentative_arg==1
								br_plotter.options.render_vertices = true;
							elseif tentative_arg==0
								br_plotter.options.render_vertices = false;
							else
								error('bad option %f for vertices',tentative_arg);
							end
						end
						
						
					case {'filename','file'}
						br_plotter.filename = command_line_options{ii+1};
						if ~ischar(br_plotter.filename)
							error('filename argument must be a filename')
						end
					case 'proj'


						tmp = command_line_options{ii+1};

						if isa(tmp,'function_handle')
							br_plotter.options.custom_projection = tmp;
							br_plotter.options.use_custom_projection = true;
						elseif strcmpi(tmp,'natural')
	
						else
							error('value for ''proj'' must be a function handle or ''natural''');
						end

					case {'colormap'}
						tmp = command_line_options{ii+1};
						if isa(tmp,'function_handle')
							br_plotter.options.colormap = tmp;
						else
							error('value for ''colormap'' must be a handle to a function generating a colormap for an integer number of colors; e.g. @jet');
						end
						
					case {'colorfn'}
						tmp = command_line_options{ii+1};
						if isa(tmp,'function_handle')
							br_plotter.options.use_colorfn = true;
							br_plotter.options.colorfn = tmp;
						else
							error('value for ''colorfn'' must be a handle to a function accepting a vector and returning a scalar');
						end
					
					
					case {'num_colors','numcolors'}
						tmp = command_line_options{ii+1};
						if ~isint(tmp)
							error('value for ''num_colors'' must be in integer');
						end
						
						br_plotter.options.num_colors = tmp;
						
					case {'mono','monocolor'}
						br_plotter.options.monocolor = true;
						
						tentative_color = command_line_options{ii+1};
						
						
						
						if ischar(tentative_color)
							switch tentative_color
								case 'r'
									br_plotter.options.monocolor_color = [1 0 0];
								case 'g'
									br_plotter.options.monocolor_color = [0 1 0];
								case 'b'
									br_plotter.options.monocolor_color = [0 0 1];
								case 'm'
									br_plotter.options.monocolor_color = [1 0 1];
								case 'c'
									br_plotter.options.monocolor_color = [0 1 1];
								case 'y'
									br_plotter.options.monocolor_color = [1 1 0];
								case 'k'
									br_plotter.options.monocolor_color = [0 0 0];
								
								otherwise
									error('input color string for mono must be one of r g b m c y k.  you can also specify a 1x3 RGB color vector');
							end
						else
							[m,n] = size(tentative_color);
							if and(m==1,n==3)
								br_plotter.options.monocolor_color = tentative_color;
							else
								error('explicit color for monocolor surfaces must be a 1x3 RGB vector');
							end
							
						end
						
						
						
						br_plotter.options.colormap = @(num_colors) repmat(br_plotter.options.monocolor_color,num_colors,1);
						
					case 'labels'
						tentative_arg = command_line_options{ii+1};
						
						if ischar(tentative_arg)
							switch tentative_arg
								case {'y','yes','true'}
									br_plotter.options.labels = true;
								case {'n','no','none','false'}
									br_plotter.options.labels = false;
								otherwise
									error('bad option %s for labels',tentative_arg);
							end
							
						else
							if tentative_arg==1
								br_plotter.options.labels = true;
							elseif tentative_arg==0
								br_plotter.options.labels = false;
							else
								error('bad option %f for labels',tentative_arg);
							end
						end
						
					case 'linestyle'
						br_plotter.options.linestyle = command_line_options{ii+1};
						br_plotter.options.use_fixed_linestyle = true;
						
						
					case 'faces'
						
						tentative_arg = command_line_options{ii+1};
						
						if ischar(tentative_arg)
							switch tentative_arg
								case {'y','yes','true'}
									br_plotter.options.render_faces = true;
								case {'n','no','none','false'}
									br_plotter.options.render_faces = false;
								otherwise
									error('bad option %s for faces',tentative_arg);
							end
							
						else
							if tentative_arg==1
								br_plotter.options.render_faces = true;
							elseif tentative_arg==0
								br_plotter.options.render_faces = false;
							else
								error('bad option %f for faces',tentative_arg);
							end
						end
						
												
					otherwise
						error('unexpected option name ''%s''',command_line_options{ii})
				end
			end
		end
		
		
		
		% uses internally set variable 'filename' to load a .mat file
		% containing data gathered previously.
		function load_data(br_plotter)
			if isempty(br_plotter.filename)
				error('unset filename in br_plotter object');
			end
			
			if isempty(dir(br_plotter.filename))
				error('nexists file with name ''%s''',br_plotter.filename);
			end
			
			file_variables = whos('-file',br_plotter.filename);
			
			if ismember('BRinfo', {file_variables.name})
				temp = load(br_plotter.filename);
				br_plotter.BRinfo = temp.BRinfo;
			else
				error('file ''%s'' does not contain variable ''BRinfo''',br_plotter.filename);
			end
			
			[br_plotter.options.containing, br_plotter.options.basename, ~] = fileparts(pwd);
			br_plotter.dimension = br_plotter.BRinfo.dimension;
		end
		
		
		function set_filename(br_plotter, new_filename)
			
			if nargin==1
				prev_filenames = dir('BRinfo*.mat');
				
				if isempty(prev_filenames)
					br_plotter.filename = uigetfile();
				else
					max_found = -1;

					for ii = 1:length(prev_filenames)
						curr_name = prev_filenames(ii).name;
						curr_num = str2double(curr_name(7:end-4));
						if max_found < curr_num
							max_found = curr_num;
						end
					end
					br_plotter.filename = ['BRinfo' num2str(max_found) '.mat'];
				end
			else	
				br_plotter.filename = new_filename;
			end
			
		end
		

		
		function load_and_render(br_plotter,source, event)
			
			[FileName,PathName,FilterIndex] = uigetfile();
			br_plotter.filename = [PathName FileName];
			
			load(br_plotter);
			plot(br_plotter);
			
		end
		
		
		function plot(br_plotter)
			
			if br_plotter.BRinfo.num_vertices==0
				warning('your decomposition contains 0 vertices.  The real part appears to be empty.')
				return;
			end
				
			setupfig(br_plotter);
			
			
			switch br_plotter.dimension
				case 1
					curve_plot(br_plotter);

				case 2	
					surface_plot(br_plotter);

				otherwise

			end

			if ~isempty(br_plotter.legend)
				render_legends(br_plotter);
			end
			
			
			br_plotter.options.plotted = 1;
			
			
			
			button_setup(br_plotter);
			
			controls(br_plotter);
			
			if br_plotter.options.autosave
				try
					save_routine(br_plotter);
				catch exception
					display(exception);
					display('saving render not completed');
				end
			end
			
			
		end
		
		
		function set_label_text_size(br_plotter,~,~,new_size)
			
			f = fieldnames(br_plotter.handles.vertex_text);
			for ii = 1:length(f)
				set(br_plotter.handles.vertex_text.(f{ii}),'FontSize',new_size);
			end
			
			
			switch br_plotter.dimension
				case 1
					
					
				case 2
					all_edge_labels = [br_plotter.handles.critcurve_labels;...
								br_plotter.handles.spherecurve_labels;...
								br_plotter.handles.midtext;...
								br_plotter.handles.crittext;...
								br_plotter.handles.singtext];
					
					set(all_edge_labels,'FontSize',new_size);
				otherwise
					
					
			end
			
			br_plotter.options.fontsizes.labels = new_size;
		end
	
		function set_legend_text_size(br_plotter,~,~,new_size)
			set(br_plotter.handles.legend,'FontSize',new_size);
			br_plotter.options.fontsizes.legend = new_size;
		end
		
		function set_axis_text_size(br_plotter,~,~,new_size)
			set(br_plotter.axes.main,'FontSize',new_size);
			br_plotter.options.fontsizes.axis = new_size;
		end
		
		
		function delete_panels(br_plotter)
			
			b = fieldnames(br_plotter.panels);
			for ii = 1:length(b)
				delete(br_plotter.panels.(b{ii}));
			end
			br_plotter.panels = [];
			
		end
		
		
		
		
		
		
		
		% declared headers for functions in other files.
		
		preprocess_data(br_plotter)

		get_indices(br_plotter)
		
		
		
		
		%functions specific to surfaces
		surface_plot(br_plotter)
		handles = plot_surface_edges(br_plotter)
		
		
		%functions specific to curves
		curve_plot(br_plotter)
		handles = plot_curve_samples(br_plotter,sampler_data,style, color)
		
		
		% common to all dimensions
		sphere_plot(br_plotter)
		plot_vertices(br_plotter)
		plot_edge_points(br_plotter)
		
		setupfig(br_plotter,varargin)
		
		create_axes(br_plotter)
		label_axes(br_plotter)
		sync_axes(br_plotter)
		
		render_legends(br_plotter)
		visibility_setup(br_plotter)
		
		
		
		% calls the initial_visibility, visibility_setup, and
		% twiddle_visibility functions
		controls(br_plotter)
		camera_setup(br_plotter)
		
		
		
		
		% setup the interactive buttons
		button_setup(br_plotter)
		set_initial_visibility(br_plotter)
		twiddle_visibility(br_plotter)
		
		
		
		
		% callback functions
		
		% for more info on associating callbacks with buttons, see e.g.
		% http://www.mathworks.com/help/matlab/matlab_oop/class-methods-for-graphics-callbacks.html
		change_alpha(br_plotter,source,event) % is a callback function
		change_text_size(br_plotter,source,event)% is a callback function
		center_camera_on_selected_point(br_plotter,source, event)
		save_routine(br_plotter,varargin) % is a callback function
		flip_switch(br_plotter,srcHandle,eventData,varargin)
		resizeui(br_plotter,srcHandle,eventData,varargin)
		
		
	
	end%re: methods
	
	

	
end
