function fv = curve_pipe_surface(BRinfo, use_sampled_data, varargin)
%  generate a pipe surface from a bertini_real decomposition of a curve
%
% the input BRinfo must be a curve.
%
% vararg-pair options - n, radius, render, write_stl
%

if BRinfo.dimension~=1
	error('the decomposition must be a curve.');
end

opt = set_options(varargin);


switch BRinfo.num_variables-1
	case 2
		variable_indices = 1:2;
	case 3
		variable_indices = 1:3;
	otherwise
		variable_indices = 1:3;
end


fv.faces = [];
fv.vertices = [];


raw_edges = BRinfo.edges;

degenerate = raw_edges(:,1)==raw_edges(:,3); %discard the degenerate edges

br_vertices_needing_convhull = [];
fv_vertices_to_convhull = {};
for edge_index = 1:BRinfo.num_edges
	
	if use_sampled_data
		indices = BRinfo.sampler_data.edge(edge_index).samples+1;
	else
		indices = BRinfo.edges(edge_index,:);
	end
	
	if indices(1)==indices(end)
		continue
	end
	
	
	edge = real([BRinfo.vertices(indices).point]);
	
	
	if BRinfo.num_variables-1==2
		edge = [edge;zeros(1,size(edge,2))]; %edge = [edge;edge(2,:)];%
	elseif BRinfo.num_variables-1>=4
		edge = edge(variable_indices,:);
	end
	

	match_left = raw_edges(:,1)==indices(1);
	match_right = raw_edges(:,3)==indices(1);
	closed_left = sum(and(match_left,~degenerate) + and(match_right,~degenerate))<=1;
	
	
	match_left = raw_edges(:,1)==indices(end);
	match_right = raw_edges(:,3)==indices(end);
	closed_right = sum(and(match_left,~degenerate) + and(match_right,~degenerate))<=1;
	
	
	if closed_left && closed_right
		closed_val = 'both';
	elseif closed_left
		closed_val = 'left';
	elseif closed_right
		closed_val = 'right';
	else
		closed_val = 'none';
	end
	
	[h,temp_fv] = pipe_surface(edge(1,:),edge(2,:),edge(3,:),'r',opt.radius,'n',opt.n,'closed',closed_val,'render',false);
	
	if sum(sum(isnan(temp_fv.vertices)))>0
		warning('pipe surface has nans');
	end
% 	edge_index
% 	fv
% 	temp_fv
% 	
% 	
	if ~closed_left
		
		%store the first n+1 vertex indices as those which need to have
		%convhull done to them.
		ind = find(br_vertices_needing_convhull==indices(1),1);
		if isempty(ind)
			br_vertices_needing_convhull(end+1) = indices(1);
			ind = length(br_vertices_needing_convhull);
			fv_vertices_to_convhull{ind} = [];
		end
		%the first n+1 vertices in a temp_fv are those for the left end's
		%circle
		fv_vertices_to_convhull{ind} = [fv_vertices_to_convhull{ind} size(fv.vertices,1)+temp_fv.left_cap];
	end
	
	if ~closed_right

		%store the first n+1 vertex indices as those which need to have
		%convhull done to them.
		ind = find(br_vertices_needing_convhull==indices(end),1);
		if isempty(ind)
			br_vertices_needing_convhull(end+1) = indices(end);
			ind = length(br_vertices_needing_convhull);
			fv_vertices_to_convhull{ind} = [];
		end
		%the first n+1 vertices in a temp_fv are those for the left end's
		%circle
		
		fv_vertices_to_convhull{ind} = [fv_vertices_to_convhull{ind} size(fv.vertices,1)+temp_fv.right_cap];
	end
	
	fv.faces = [fv.faces;temp_fv.faces+size(fv.vertices,1)];
	fv.vertices = [fv.vertices;temp_fv.vertices];
	
end


convhulls.vertices = fv.vertices;
convhulls.faces = [];

for ii = 1:length(fv_vertices_to_convhull)
	local_indices = fv_vertices_to_convhull{ii};
	
	try
		q = convhulln(fv.vertices(local_indices,:)); %(end:-1:1)
		q = local_indices(q);
		q = q(:,[3 2 1]);
		fv.faces = [fv.faces;q];
		convhulls.faces = [convhulls.faces;q];
	catch
%		warning('computation of convex hull for index %i failed',ii);
	end
	
	
end



if opt.render
	h = patch(convhulls);
	set(h,'FaceAlpha',0.2);
	set(gca,'DataAspectRatio',[1 1 1]);
	set(h,'FaceColor',[0 1 0]);
	set(h,'EdgeColor',[0 0.5 0]);
	cameratoolbar
	rotate3d on 
end

if opt.write_to_stl
	fv_to_stl(fv)
end


end





function opt = set_options(command_line_options)

if mod(length(command_line_options),2)~=0
	error('must have option-value pairs');
end

opt.radius = 0.15;
opt.n = 31;
opt.render = true;
opt.write_to_stl = false;

for ii = 1:2:length(command_line_options)-1
	val = command_line_options{ii+1};
	
	switch command_line_options{ii}
		case 'radius'
			opt.radius = val;
		case 'n'
			opt.n = val;
		case 'render'
			opt.render = val;
		case 'write_stl'
			opt.write_to_stl = val;
		otherwise
			error('bad option %s',command_line_options{ii});
	end
end

end




