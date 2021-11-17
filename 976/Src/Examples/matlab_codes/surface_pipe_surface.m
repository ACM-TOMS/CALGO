function fv = surface_pipe_surface(BRinfo, varargin)
% generate a pipe surface from a bertini_real decomposition of a surface
%
%
%


if BRinfo.dimension~=2
	error('input decomposition must be a surface');
end

opt = set_options(varargin, BRinfo);



fv.faces = [];
fv.vertices = [];


%get only the relevant coordinates in the points in the vertex set.
pc_vertices = BRinfo.vertices;
for ii = 1:length(pc_vertices)
	pc_vertices(ii).point = pc_vertices(ii).point(1:BRinfo.num_variables-1);
end



if opt.crit
	curve = BRinfo.crit_curve;
	curve.vertices = pc_vertices;
	
	temp_fv = curve_pipe_surface(curve,opt.use_sampled,...
			'radius',opt.radius,'n',opt.n,'render',opt.render_curves,'write_stl',false);
	
	fv.faces = [fv.faces; temp_fv.faces+size(fv.vertices,1)];
	fv.vertices = [fv.vertices; temp_fv.vertices];
end

if opt.render_curves
	hold on
end

if opt.sphere
	curve = BRinfo.sphere_curve;
	curve.vertices = pc_vertices;
	
	temp_fv = curve_pipe_surface(curve,opt.use_sampled,...
			'radius',opt.radius,'n',opt.n,'render',opt.render_curves,'write_stl',false);
	
	fv.faces = [fv.faces; temp_fv.faces+size(fv.vertices,1)];
	fv.vertices = [fv.vertices; temp_fv.vertices];
end

if opt.sing
	for ii = 1:length(BRinfo.singular_curves)
		curve = BRinfo.singular_curves{ii};
		curve.vertices = pc_vertices;

		temp_fv = curve_pipe_surface(curve,opt.use_sampled,...
			'radius',opt.radius,'n',opt.n,'render',opt.render_curves,'write_stl',false);

		fv.faces = [fv.faces; temp_fv.faces+size(fv.vertices,1)];
		fv.vertices = [fv.vertices; temp_fv.vertices];
	end
end

if opt.midslice
	for ii = 1:length(BRinfo.midpoint_slices)
		curve = BRinfo.midpoint_slices{ii};
		curve.vertices = pc_vertices;
		temp_fv = curve_pipe_surface(curve,opt.use_sampled,...
			'radius',opt.radius,'n',opt.n,'render',opt.render_curves,'write_stl',false);

		fv.faces = [fv.faces; temp_fv.faces+size(fv.vertices,1)];
		fv.vertices = [fv.vertices; temp_fv.vertices];
	end
end

if opt.critslice
	for ii = 1:length(BRinfo.critpoint_slices)
		curve = BRinfo.critpoint_slices{ii};
		curve.vertices = pc_vertices;
		temp_fv = curve_pipe_surface(curve,opt.use_sampled,...
			'radius',opt.radius,'n',opt.n,'render',opt.render_curves,'write_stl',false);

		fv.faces = [fv.faces; temp_fv.faces+size(fv.vertices,1)];
		fv.vertices = [fv.vertices; temp_fv.vertices];
	end
end

if opt.render
	h = patch(fv);
	set(h,'FaceAlpha',0.1);
	set(h,'EdgeAlpha',0.1);
	set(h,'FaceColor',[0 0 1]);
	set(h,'EdgeColor',[0 0 0.5]);
	set(gca,'DataAspectRatio',[1 1 1])
	rotate3d on 
end

hold off

end


function opt = set_options(command_line_options, BRinfo)

if mod(length(command_line_options),2)~=0
	error('must have option-value pairs');
end

opt.crit = true;
opt.sing = true;
opt.sphere = true;
opt.midslice = true;
opt.critslice = true;

opt.radius = 0.15;
opt.n = 31;

opt.render_curves = true;
opt.render = true;
opt.write_to_stl = false;

if ~isempty(BRinfo.sampler_data)
	opt.use_sampled = true;
else
	opt.use_sampled = false;
end

for ii = 1:2:length(command_line_options)-1
	val = command_line_options{ii+1};
	
	switch command_line_options{ii}
		case 'crit'
			opt.crit = val;
		case 'sing'
			opt.sing = val;
		case 'sphere'
			opt.sphere = val;
		case 'midslice'
			opt.midslice = val;
		case 'critslice'
			opt.critslice = val;
		case 'n'
			opt.n = val;
		case 'radius'
			opt.radius = val;
		case 'render'
			opt.render = val;
		case 'write_stl'
			opt.write_to_stl = val;
		case 'sampled'
			opt.use_sampled = val;
		otherwise
			error('bad option %s',command_line_options{ii});
	end
end

end
