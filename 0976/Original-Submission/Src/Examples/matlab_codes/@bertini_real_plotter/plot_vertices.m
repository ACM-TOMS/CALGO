function br_plotter = plot_vertices(br_plotter)

markersize = br_plotter.options.markersize;
ind = br_plotter.indices;
curr_axis = br_plotter.axes.main;
fontsize = br_plotter.options.fontsizes.labels;

%these should be made programmatic, after modding bertini_real to generate
%the table of them for us.

[names, types] = vertex_type_indices(br_plotter);

markers = {'x','o','s','d','^','v','>','<','p','h'};
colors = [1 0 0;lines(length(names)-1)];

for ii = 1:length(names)
	indexed_markers.(names{ii}).marker = markers{ii};
	indexed_markers.(names{ii}).color = colors(ii,:);
end







if br_plotter.options.labels
	labels = cell(br_plotter.BRinfo.num_vertices,1);
end

br_plotter.fv.vertices = zeros(br_plotter.BRinfo.num_vertices,length(ind));



unpacked_vertex_types = zeros(br_plotter.BRinfo.num_vertices, 1);
for ii=1:br_plotter.BRinfo.num_vertices
	
    br_plotter.fv.vertices(ii,:) = transpose(real(br_plotter.BRinfo.vertices(ii).point(ind)));
	
	unpacked_vertex_types(ii) = br_plotter.BRinfo.vertices(ii).type;
	
	if br_plotter.options.labels
		labels{ii} = [ '    ' num2str(ii-1)];
	end
end

using_bitor = isfield(br_plotter.BRinfo,'vertex_types');
if using_bitor
	has_type = @(x,y) bitand(x,y)>0;
else
	has_type = @(x,y) x==y;
end

if br_plotter.options.render_vertices
	rendered_counter = 0;
	for ii = 1:length(types)

		curr_name = names{ii};
		
		curr_indices_logical = has_type(unpacked_vertex_types, types(ii));

		if sum(curr_indices_logical)==0
			continue;
		else
			rendered_counter = rendered_counter+1;
		end

		curr_points = br_plotter.fv.vertices(curr_indices_logical,:);


		h = main_plot_function(curr_points,   1:length(ind), curr_axis);


		set(h,'LineStyle','none');
		set(h,'Marker', indexed_markers.(curr_name).marker,'MarkerSize',markersize,'MarkerFaceColor',0.6*indexed_markers.(curr_name).color,'MarkerEdgeColor',indexed_markers.(curr_name).color);
		set(h,'Color',colors(ii,:));

		local_catname = sprintf('ind_%d',types(ii));



		br_plotter.handles.vertices.(names{ii}) = h;



		br_plotter.legend.vertices.handles(rendered_counter) = h;
		br_plotter.legend.vertices.types{rendered_counter} = names{ii};
		br_plotter.legend.vertices.text{rendered_counter} = ['point type ' names{ii}];


		if br_plotter.options.labels
			switch length(ind)

				case 2
					br_plotter.handles.vertex_text.(names{ii}) = text(curr_points(:,1), curr_points(:,2), labels(curr_indices_logical),...
						'HorizontalAlignment','left','VerticalAlignment','bottom',...
						'FontSize',fontsize,'Parent',curr_axis,'Color',colors(ii,:));
				case 3
					br_plotter.handles.vertex_text.(names{ii}) = text(curr_points(:,1), curr_points(:,2), curr_points(:,3), labels(curr_indices_logical),...
						'HorizontalAlignment','left','VerticalAlignment','bottom',...
						'FontSize',fontsize,'Parent',curr_axis,'Color',colors(ii,:));
				otherwise
			end
		end

	end
end

end% re function


function [names, types] = vertex_type_indices(br_plotter)


if br_plotter.BRinfo.run_metadata.version.number < 104
	
	names = {'UNSET', 'CRITICAL', 'SEMICRITICAL', ...
			'MIDPOINT', 'ISOLATED', 'NEW', ...
			'CURVE_SAMPLE_POINT', 'SURFACE_SAMPLE_POINT', 'REMOVED', ...
			'PROBLEMATIC'};
	catnames = cell(length(names),1);
	curr_index = 100; %sadly set manually...  this corresponds to header files in bertini_real.  

	for ii = 1:length(names)
		catnames{ii} = ['ind_' num2str(curr_index)];
		available_types.(catnames{ii}) = names{ii};
		curr_index = curr_index+1;
	end
	types = 100:109; %initialize to blank  eww i hate these magic constants.  sorry. the code in 1.4 addresses this stupid problem
else
	names = br_plotter.BRinfo.vertex_types.names;
	types = br_plotter.BRinfo.vertex_types.nums;
end

end
