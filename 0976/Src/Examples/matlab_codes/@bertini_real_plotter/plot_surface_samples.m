

function plot_surface_samples(br_plotter)


if br_plotter.options.monocolor
	sampler_faces = plot_surf_samples_monocolor(br_plotter);
elseif br_plotter.options.use_colorfn
	sampler_faces = plot_surf_samples_colorfn(br_plotter);
else
	sampler_faces = plot_surf_samples_multicolor(br_plotter);
end


if ~isempty(sampler_faces)
	br_plotter.fv.faces = sampler_faces;
end

end



function sampler_faces = plot_surf_samples_colorfn(br_plotter)

total_num_faces = 0;
for ii = 1:length(br_plotter.BRinfo.sampler_data)
	total_num_faces = total_num_faces+size(br_plotter.BRinfo.sampler_data{ii},1);
end

sampler_faces = zeros(total_num_faces,3);




ind_so_far = 1;

tmp_fv = br_plotter.fv;

sample_cdata = zeros(size(tmp_fv.vertices,1),1);
for ii = 1:size(tmp_fv.vertices,1)
	sample_cdata(ii) = br_plotter.options.colorfn(tmp_fv.vertices(ii,:));
end


if ~isempty(br_plotter.BRinfo.sampler_data)
	for ii = 1:length(br_plotter.BRinfo.sampler_data)
		
		
		
		tmp_fv.faces = br_plotter.BRinfo.sampler_data{ii}+1;  %shift by 1 here to correct for 1-based indexing
		tmp_fv.faces(any(br_plotter.fv.faces<=0,2),:) = []; % omit problematic faces.
		
		num_this_time = size(tmp_fv.faces,1);
		
		%plot the BR face
		h = patch(tmp_fv);
		set(h,...
			'FaceVertexCData', sample_cdata,...
			'FaceColor', 'interp',...
			'FaceAlpha',br_plotter.options.face_alpha,'EdgeColor','none','EdgeAlpha',br_plotter.options.edge_alpha);
		br_plotter.handles.surface_samples(ii) = h;
		
		% add the faces to the total face blabal
		sampler_faces(ind_so_far:ind_so_far+num_this_time-1,:) = tmp_fv.faces;
		ind_so_far = ind_so_far + num_this_time;
	end
end


sampler_faces = sampler_faces(1:ind_so_far-1,:);


end



function sampler_faces = plot_surf_samples_multicolor(br_plotter)

colors = br_plotter.options.colormap(length(br_plotter.BRinfo.sampler_data));

total_num_faces = 0;
for ii = 1:length(br_plotter.BRinfo.sampler_data)
	total_num_faces = total_num_faces+size(br_plotter.BRinfo.sampler_data{ii},1);
end

sampler_faces = zeros(total_num_faces,3);




ind_so_far = 1;

tmp_fv = br_plotter.fv;

if ~isempty(br_plotter.BRinfo.sampler_data)
	for ii = 1:length(br_plotter.BRinfo.sampler_data)
		
		
		
		tmp_fv.faces = br_plotter.BRinfo.sampler_data{ii}+1;
		tmp_fv.faces(any(br_plotter.fv.faces<=0,2),:) = []; % omit problematic faces.
		
		num_this_time = size(tmp_fv.faces,1);
		
		%plot the BR face
		h = patch(tmp_fv);
		set(h,'FaceColor',colors(ii,:),'FaceAlpha',br_plotter.options.face_alpha,'EdgeColor',0.8*colors(ii,:),'EdgeAlpha',br_plotter.options.edge_alpha);
		br_plotter.handles.surface_samples(ii) = h;
		
		% add the faces to the total face blabal
		sampler_faces(ind_so_far:ind_so_far+num_this_time-1,:) = tmp_fv.faces;
		ind_so_far = ind_so_far + num_this_time;
	end
end


sampler_faces = sampler_faces(1:ind_so_far-1,:);


end





function sampler_faces = plot_surf_samples_monocolor(br_plotter)



total_num_faces = 0;
for ii = 1:length(br_plotter.BRinfo.sampler_data)
	total_num_faces = total_num_faces+size(br_plotter.BRinfo.sampler_data{ii},1);
end

sampler_faces = zeros(total_num_faces,3);


br_plotter.handles.surface_samples = [];

ind_so_far = 1;

tmp_fv = br_plotter.fv;

if ~isempty(br_plotter.BRinfo.sampler_data)
	for ii = 1:length(br_plotter.BRinfo.sampler_data)
		
		
		
		tmp_fv.faces = br_plotter.BRinfo.sampler_data{ii}+1;
		tmp_fv.faces(any(br_plotter.fv.faces<=0,2),:) = []; % omit problematic faces.
		
		num_this_time = size(tmp_fv.faces,1);
		
		% add the faces to the total face blabal
		sampler_faces(ind_so_far:ind_so_far+num_this_time-1,:) = tmp_fv.faces;
		ind_so_far = ind_so_far + num_this_time;
	end
end



sampler_faces = sampler_faces(1:ind_so_far-1,:);

tmp_fv.faces = sampler_faces;

%plot the BR samples in a single color
thecolor = br_plotter.options.monocolor_color;

h = patch(tmp_fv);
set(h,'FaceColor',thecolor,'FaceAlpha',br_plotter.options.face_alpha,'EdgeColor',thecolor,'EdgeAlpha',br_plotter.options.edge_alpha);
br_plotter.handles.surface_samples(1) = h;
		
		

end
