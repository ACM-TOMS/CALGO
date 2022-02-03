

function br_plotter = plot_surface_edges(br_plotter)




handle_counter = 0;


if br_plotter.BRinfo.crit_curve.num_edges>0
	
	if isfield(br_plotter.options,'use_fixed_linestyle')
			style = br_plotter.options.linestyle;
		else
			style = '-';
	end
		
	[br_plotter.handles.critcurve, ...
		br_plotter.handles.refinements.critcurve, ...
		br_plotter.handles.critcurve_labels] ...
		  = plot_subcurve(br_plotter,br_plotter.BRinfo.crit_curve,'critcurve',style,'r');
	
	

	if ~isempty(br_plotter.handles.critcurve)
		handle_counter = handle_counter+1;
		br_plotter.legend.surface_edges.handles(handle_counter) = br_plotter.handles.critcurve(1);
		br_plotter.legend.surface_edges.text{handle_counter} = 'crit curve';
	end
end










if isfield(br_plotter.BRinfo,'sphere_curve')
	if br_plotter.BRinfo.sphere_curve.num_edges>0
		
		if isfield(br_plotter.options,'use_fixed_linestyle')
			style = br_plotter.options.linestyle;
		else
			style = '-.';
		end
		
		[br_plotter.handles.spherecurve, ...
			br_plotter.handles.refinements.spherecurve, ...
			br_plotter.handles.spherecurve_labels] ...
			  = plot_subcurve(br_plotter,br_plotter.BRinfo.sphere_curve,'spherecurve',style,'c');
	
		if ~isempty(br_plotter.handles.spherecurve)
			handle_counter = handle_counter+1;
			br_plotter.legend.surface_edges.handles(handle_counter) = br_plotter.handles.spherecurve(1);
			br_plotter.legend.surface_edges.text{handle_counter} = 'sphere curve';
		end
	end
end



num_midslices = length(br_plotter.BRinfo.midpoint_slices);
num_crit_slices = length(br_plotter.BRinfo.critpoint_slices);

colors = br_plotter.options.colormap(2*num_midslices+1);






firstone = 1;

if isfield(br_plotter.options,'use_fixed_linestyle')
	style = br_plotter.options.linestyle;
else
	style = ':';
end
		
for kk = 1:num_crit_slices
	color_index = 1+2*(kk-1);
	
	[handie_mc_handhand, refinement_handles, label_handles] ...
		= plot_subcurve(br_plotter,br_plotter.BRinfo.critpoint_slices{kk},sprintf('crit.%d.',kk-1),style,'m',colors(color_index,:));
	
	br_plotter.handles.crittext = [br_plotter.handles.crittext;label_handles];
	br_plotter.handles.critslices = [br_plotter.handles.critslices;handie_mc_handhand];
	br_plotter.handles.refinements.critslice = [br_plotter.handles.refinements.critslice;refinement_handles];
	
	if firstone
		if ~isempty(handie_mc_handhand)
			handle_counter = handle_counter+1;
			br_plotter.legend.surface_edges.handles(handle_counter) = handie_mc_handhand(1);
			br_plotter.legend.surface_edges.text{handle_counter} = 'critslices';
			firstone = 0;
		end
	end
	
	
end






if isfield(br_plotter.options,'use_fixed_linestyle')
	style = br_plotter.options.linestyle;
else
	style = '--';
end
		

added = false;
for kk = 1:length(br_plotter.BRinfo.midpoint_slices)
	if strcmp( br_plotter.BRinfo.midpoint_slices{kk}.inputfilename,'unset')
		continue;
	end
	color_index = 2*(kk);
	[handie_mc_handhand, refinement_handles, label_handles] ...
		= plot_subcurve(br_plotter,br_plotter.BRinfo.midpoint_slices{kk},sprintf('mid.%d.',kk-1),style,'g',colors(color_index,:));
	
	br_plotter.handles.midtext = [br_plotter.handles.midtext;label_handles];
	br_plotter.handles.midslices = [br_plotter.handles.midslices;handie_mc_handhand];
	br_plotter.handles.refinements.midslice = [br_plotter.handles.refinements.midslice;refinement_handles];
	
	if and(added==false,~isempty(handie_mc_handhand))
		handle_counter = handle_counter+1;
		br_plotter.legend.surface_edges.handles(handle_counter) = handie_mc_handhand(1);
		br_plotter.legend.surface_edges.text{handle_counter} = 'midslices';
		added = true;
	end
	
	
end









if isfield(br_plotter.options,'use_fixed_linestyle')
	style = br_plotter.options.linestyle;
else
	style = ':';
end
		
if isfield(br_plotter.BRinfo,'singular_curves')
	
	
	firstone = 1;

	for kk = 1:length(br_plotter.BRinfo.singular_curves)
		
		
		[handie_mc_handhand, refinement_handles, label_handles] = plot_subcurve(br_plotter,br_plotter.BRinfo.singular_curves{kk},sprintf('sing.%d.',kk),style,'b');
	

		br_plotter.handles.singtext = [br_plotter.handles.singtext;label_handles];
		br_plotter.handles.singular_curves = [br_plotter.handles.singular_curves;handie_mc_handhand];
		br_plotter.handles.refinements.singularcurve = [br_plotter.handles.refinements.singularcurve;refinement_handles];

		if firstone
			if ~isempty(handie_mc_handhand)
				handle_counter = handle_counter+1;
				br_plotter.legend.surface_edges.handles(handle_counter) = handie_mc_handhand(1);
				br_plotter.legend.surface_edges.text{handle_counter} = 'singular_curves';
				firstone = 0;
			end
		end
	end
end




end









function [edge_handles, refinement_handles, text_handle] = plot_subcurve(br_plotter,curve,name,style,text_color,desiredcolor)


curr_axes = br_plotter.axes.main;


edge_handles = [];
refinement_handles = [];
text_handle = [];



num_nondegen = 0;
nondegen_edge_indices = zeros(curve.num_edges,1);
for ii = 1:curve.num_edges
	if curve.edges(ii,1)~=curve.edges(ii,3)
		num_nondegen = num_nondegen+1;
		nondegen_edge_indices(num_nondegen) = ii;
	end
end
nondegen_edge_indices = nondegen_edge_indices(1:num_nondegen);



if num_nondegen==0
	return;
end



if nargin==6
	colors = repmat(desiredcolor,[num_nondegen 1]);
else
	colors = 0.8*br_plotter.options.colormap(num_nondegen);
end


edge_handles = zeros(num_nondegen,1);

if br_plotter.options.labels
	text_locations = zeros(num_nondegen,3);
	text_labels = cell(num_nondegen,1);
end


for ii =1:num_nondegen
	curr_edge_index = nondegen_edge_indices(ii);
	curr_edge = curve.edges(curr_edge_index,:);
	

	curve_edge_points = br_plotter.fv.vertices(curr_edge,:);

	h = plot3(curve_edge_points(:,1),curve_edge_points(:,2),curve_edge_points(:,3),'Parent',curr_axes);

	set(h,'Color',colors(ii,:));
	set(h,'LineStyle',style,'LineWidth',br_plotter.options.line_thickness);


	
	edge_handles(ii) = h;
	
	if br_plotter.options.labels
		text_locations(ii,:) = curve_edge_points(2,:);
		text_labels{ii} = sprintf('%s %d  ',name,curr_edge_index-1);
	end
	
end

if br_plotter.options.labels
	text_handle = text(text_locations(:,1),text_locations(:,2),text_locations(:,3), text_labels,...
				'HorizontalAlignment','right',...
				'FontSize',br_plotter.options.fontsizes.labels,...
				'Parent',curr_axes,'Color',text_color);
end


if ~isempty(br_plotter.BRinfo.sampler_data)
	if nargin==6
		refinement_handles = plot_curve_samples(br_plotter,curve.sampler_data,style, desiredcolor);
	else
		refinement_handles = plot_curve_samples(br_plotter,curve.sampler_data,style);
	end
end


		

end
