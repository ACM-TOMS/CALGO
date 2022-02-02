function gather_br_samples()



[dirname,dimension] = parse_dirname;

locations = strfind(dirname,'/');
if length(locations)>0
	dirname = dirname(locations(end)+1:end);
end

BRinfo.dirname = dirname;
BRinfo.dimension = dimension;







BRinfo = gather_vertices(BRinfo);



switch BRinfo.dimension
	
	case 1
		
		[BRinfo] = gather_curve(BRinfo.dirname, BRinfo);

	case 2
		
		[BRinfo] = gather_surface(BRinfo.dirname, BRinfo);
		
	otherwise
		
		
	
end

BRinfo.run_metadata = gather_run_metadata(BRinfo.dirname);
BRinfo.vertex_types = gather_vertex_types(BRinfo.dirname);

tmpnames = get_names(BRinfo.num_variables);
BRinfo.var_names = tmpnames(2:end);

display('done gathering decomposition.');

filename = generate_filename();

display('writing data to file.');
save(filename,'BRinfo');

display(['file saved to ' filename]);
end%re: function


function filename = generate_filename()

prev_filenames = dir('BRinfo*.mat');
max_found = -1;

for ii = 1:length(prev_filenames)
	curr_name = prev_filenames(ii).name;
	curr_num = str2num(curr_name(7:end-4));
	if max_found < curr_num
		max_found = curr_num;
	end
	
end

filename = ['BRinfo' num2str(max_found+1) '.mat'];
end


function md = gather_run_metadata(dirname)
	fname = [dirname '/run_metadata'];
	if ~isempty(dir(fname))

		fin = fopen(fname,'r');

		md.version.string = fgetl(fin);
		md.dir.string = fgetl(fin);
		md.time.string = fgetl(fin);
		md.parallel.numprocs = fscanf(fin,'%i');
		fclose(fin);

		pds = strfind(md.version.string, '.');

		md.version.major = str2num(md.version.string(1:pds(1)));
		md.version.minor = str2num(md.version.string(pds(1)+1:pds(2)));
		md.version.subminor = str2num(md.version.string(pds(2)+1:end));
		
		md.version.number = 100*md.version.major + md.version.minor + 0.01 * md.version.subminor;
		
	else
		
		md.version.string = 'earlier than 1.4';
		md.version.major = 1;
		md.version.minor = 3;
		md.version.subminor = 0;
		md.version.number = 130;
	end
end


function md = gather_vertex_types(dirname)
	fname = [dirname '/vertex_types'];
	if ~isempty(dir(fname))

		fin = fopen(fname,'r');
		md.num_types = fscanf(fin, '%i\n',[1 1]);
		md.names = cell(md.num_types,1);
		md.nums = zeros(md.num_types,1);
		for ii = 1:md.num_types
			type = fscanf(fin,'%s ',[1 1]);
			num = fscanf(fin,'%i');
			
			md.by_type.(type) = num;
			md.names{ii} = type;
			md.nums(ii) = num;
		end
		fclose(fin);
		
	else
		error('no ''vertex_types'' file found.  please use the version of the matlab code which came with your bertini real version');
	end
	
end
function input = read_input(dirname,decomp)
filename = [dirname '/' decomp.inputfilename];
if ~isempty(dir(filename))
	input = fileread(filename);
else
	display(sprintf('did not find input file at %s',dirname));
	input = '';
end


end



function [BRinfo] = gather_surface(dirname, BRinfo)

BRinfo = parse_decomp(dirname,BRinfo);


fin = fopen([dirname '/S.surf'],'r');

[BRinfo.num_faces] = fscanf(fin,'%i',[1 1]);
[BRinfo.num_edges] = fscanf(fin,'%i',[1 1]);
[BRinfo.num_midpoint_slices] = fscanf(fin,'%i',[1 1]);
[BRinfo.num_critpoint_slices] = fscanf(fin,'%i',[1 1]);
BRinfo.num_singular_curves = fscanf(fin,'%i',[1 1]);
BRinfo.singular_curve_multiplicities = [];
for ii = 1:BRinfo.num_singular_curves
    BRinfo.singular_curve_multiplicities(ii,1:2) = fscanf(fin,'%i',[1 2]);
end


fclose(fin);



[BRinfo.faces] = gather_faces(dirname);

if BRinfo.num_midpoint_slices>0
	for ii =1:BRinfo.num_midpoint_slices
		a = gather_curve([dirname '/curve_midslice_' num2str(ii-1) ], []);
		BRinfo.midpoint_slices{ii} = a;
	end
else
	BRinfo.midpoint_slices = [];
end

if BRinfo.num_critpoint_slices>0
	for ii =1:BRinfo.num_critpoint_slices
		a = gather_curve([dirname '/curve_critslice_' num2str(ii-1) ],[]);
		[BRinfo.critpoint_slices{ii}] = a;
	end
else
	BRinfo.critpoint_slices = [];
end
a = gather_curve([dirname '/curve_crit'],[]);
[BRinfo.crit_curve] = a;


[BRinfo.sphere_curve] = gather_curve([dirname '/curve_sphere'],[]);

if BRinfo.num_singular_curves>0

	for ii = 1:BRinfo.num_singular_curves
		a = gather_curve([dirname '/curve_singular_mult_' num2str(BRinfo.singular_curve_multiplicities(ii,1)) '_' num2str(BRinfo.singular_curve_multiplicities(ii,2))],[]);
		
		BRinfo.singular_curves{ii} = a;
		BRinfo.singular_names{ii} = a.inputfilename;
	end
else
	BRinfo.singular_curves = [];
end

BRinfo.sampler_data = gather_surface_samples(dirname);


BRinfo.input = read_input(dirname,BRinfo);
end



function [faces] = gather_faces(dirname)

faces = [];
fid = fopen([dirname '/F.faces'],'r');

num_faces = fscanf(fid,'%i\n',[1 1]);

for ii = 1:num_faces
	[faces(ii).midpoint] = fscanf(fid,'%i', [1 1]);
	faces(ii).midslice_index = fscanf(fid,'%i\n', [1 1]);
	
	faces(ii).top = fscanf(fid,'%i',[1 1]);
	faces(ii).bottom = fscanf(fid,'%i\n',[1 1]);
	
	faces(ii).system_top = fscanf(fid,'%s',[1 1]);
	faces(ii).system_bottom = fscanf(fid,'%s\n',[1 1]);
	
	
	faces(ii).num_left = fscanf(fid,'%i\n',[1 1]); 
	faces(ii).left = zeros(1,faces(ii).num_left); % preallocate
	for jj = 1:faces(ii).num_left
		faces(ii).left(jj) = fscanf(fid,'%i\n',[1 1]);
	end
	
	[faces(ii).num_right] = fscanf(fid,'%i\n',[1 1]);
	faces(ii).right = zeros(1,faces(ii).num_right);
	for jj = 1:faces(ii).num_right
		faces(ii).right(jj) = fscanf(fid,'%i\n',[1 1]);
	end
	

end

fclose(fid);

end


function BRinfo = parse_decomp(dirname,BRinfo)


if isempty(dir([dirname '/' 'decomp']))
	display(sprintf('did not find decomp at %s',dirname));
	BRinfo = {};
	return
end

fid = fopen([dirname '/' 'decomp'],'r');

BRinfo.inputfilename = fscanf(fid,'%s\n',[1 1]);


BRinfo.num_variables = fscanf(fid,'%i',[1 1]);
BRinfo.dimension = fscanf(fid,'%i',[1 1]);




BRinfo.pi = zeros(BRinfo.num_variables,BRinfo.dimension);

for jj = 1:BRinfo.dimension
    tmp1 = fscanf(fid,'%i',[1 1]);
    tmp_num_vars = tmp1;
for ii = 1:tmp_num_vars
	tmp1 = fscanf(fid,'%e',[1 1]);
	tmp2 = fscanf(fid,'%e\n',[1 1]); % this should scan right over the empty line separating them
	BRinfo.pi(ii,jj) = tmp1+1i*tmp2;
end
end


BRinfo.pi=BRinfo.pi(2:end,:);


BRinfo.num_patches = fscanf(fid,'%i',[1 1]);
for jj = 1:BRinfo.num_patches
    BRinfo.patch.sizes(jj) = fscanf(fid,'%i',[1 1]);
    for ii = 1:BRinfo.patch.sizes(jj)
        tmp1 = fscanf(fid,'%e',[1 1]);
        tmp2 = fscanf(fid,'%e\n',[1 1]); % this should scan right over the empty line separating them
        BRinfo.patch.vectors{jj}(ii) = tmp1 + 1i*tmp2;
    end
end

tmp1 = fscanf(fid,'%e',[1 1]);
tmp2 = fscanf(fid,'%e\n',[1 1]); % this should scan right over the empty line separating them
BRinfo.radius = tmp1+1i*tmp2;
    
center_size = fscanf(fid,'%i',[1 1]);
BRinfo.center = zeros(1,center_size);
for ii = 1:center_size
    tmp1 = fscanf(fid,'%e',[1 1]);
    tmp2 = fscanf(fid,'%e\n',[1 1]); % this should scan right over the empty line separating them
    BRinfo.center(ii) = tmp1 + 1i*tmp2;
end
fclose(fid);





end


function [dirname,dimension] = parse_dirname()


if isempty('Dir_Name')
	error('no file ''Dir_Name''.  please run bertini_real');
end

fid = fopen('Dir_Name','r');
dirname = fscanf(fid,'%s',[1 1]);
MPtype = fscanf(fid,'%i',[1 1]);
dimension = fscanf(fid,'%i',[1 1]);

fclose(fid);


end


function [sampler_data] = gather_surface_samples(dirname)
sampler_data = [];

if isempty(dir([dirname '/samp.surfsamp']))
    return;
end

fid = fopen([dirname '/samp.surfsamp']);

num_faces = fscanf(fid,'%i',[1 1]);
sampler_data = cell(num_faces,1);

for ii = 1:num_faces
    
    
    num_triangles_this_edge = fscanf(fid,'%i',[1 1]);
    
    sampler_data{ii} = zeros(num_triangles_this_edge,3);
    for jj = 1:num_triangles_this_edge
        temp_triangle = fscanf(fid,'%i',[1 3]);
        sampler_data{ii}(jj,:) = temp_triangle;
    end
end


fclose(fid);



end

function [curve] = gather_curve(dirname, curve)

curve = parse_decomp(dirname,curve);

filename = [dirname '/E.edge'];

if isempty(dir(filename))
	display(sprintf('did not find edge file at %s',dirname));
	curve.num_edges = 0;
	curve.edges = [];
else
	fid = fopen(filename,'r');

	curve.num_edges = fscanf(fid,'%i\n',[1 1]);


	curve.edges = zeros(curve.num_edges,3);  % 3, as left, mid, right
	for ii = 1:curve.num_edges
		tmp = fscanf(fid,'%i',[1 3]);
		curve.edges(ii,:) = tmp+1;
	end

	fclose(fid);
	
end




filename = [dirname '/' 'samp.curvesamp'];


if isempty(dir(filename))
	curve.sampler_data = [];
else
	fid = fopen(filename,'r');
	curve.num_edges = fscanf(fid,'%i\n\n',[1 1]);
	curve.sampler_data.sample_sizes = zeros(curve.num_edges,1);
    curve.sampler_data.edge = [];
	for ii = 1:curve.num_edges
		num_samples = fscanf(fid,'%i\n',[1 1]);
		curve.sampler_data.sample_sizes(ii) = num_samples;
		for jj=1:num_samples
			tmp_ind = fscanf(fid,'%i\n',[1 1]);
			curve.sampler_data.edge(ii).samples(jj) = tmp_ind;
		end
	end
	fclose(fid);
end


filename = [dirname '/' 'curve.cnums'];
if isempty(dir(filename))
	curve.cycle_numbers = [];
else
    curve = gather_curve_cycle_numbers(filename,curve);
end


curve.input = read_input(dirname,curve);
end

function curve = gather_curve_cycle_numbers(filename, curve)

    fid = fopen(filename,'r');
    
    tmp_num_edges = fscanf(fid,'%i',[1 1]);
    if tmp_num_edges ~= curve.num_edges
        error('curve cycle numbers does not have same number of entries as the curve has edges (%i!=%i)',tmp_num_edges,curve.num_edges);
    end
    
    curve.cycle_numbers = zeros(curve.num_edges,2);  % 2, as left, right
    for ii = 1:curve.num_edges
        tmp = fscanf(fid,'%i',[1 2]);
        curve.cycle_numbers(ii,:) = tmp;
    end
    fclose(fid);
    
end


function [BRinfo] = gather_vertices(BRinfo)

fname = 'V.vertex';
if ~isempty(dir(sprintf('%s/V_samp.vertex',BRinfo.dirname)))
    fname = 'V_samp.vertex';
end

fid = fopen(sprintf('%s/%s',BRinfo.dirname,fname));
BRinfo.num_vertices = fscanf(fid,'%i',[1 1]);
num_projections = fscanf(fid,'%i',[1 1]);
num_natural_vars = fscanf(fid,'%i',[1 1]);
num_filenames = fscanf(fid,'%i',[1 1]);
BRinfo.vertices = repmat(struct('point',[]),[1 BRinfo.num_vertices]);



for ii = 1:num_projections
	for jj = 1:num_natural_vars
		tmp = fscanf(fid,'%e %e\n',[1 2]);
	end
end

BRinfo.input_filenames = cell(1,num_filenames);
for ii = 1:num_filenames
    nchars = fscanf(fid,'%i\n',[1 1]);
   BRinfo.input_filenames{ii} = fscanf(fid,'%s\n',[1 1]); 
end

for ii = 1:BRinfo.num_vertices
	tmp_num_variables = fscanf(fid,'%i',[1 1]); % number variables
	tmpvertex = zeros(tmp_num_variables,1);
	for jj = 1:tmp_num_variables
		tmp = fscanf(fid,'%e %e\n',[1 2]);
		tmpvertex(jj) = tmp(1)+1i*tmp(2);
	end
	
	
	BRinfo.vertices(ii).point = [dehomogenize(tmpvertex(1:num_natural_vars));tmpvertex(num_natural_vars+1:end)];
	
	num_proj_vals = fscanf(fid,'%i\n',[1 1]);
	for jj = 1:num_proj_vals
		tmp = fscanf(fid,'%e %e\n',[1 2]);
		BRinfo.vertices(ii).projection_value(jj) = tmp(1)+1i*tmp(2);
	end
    BRinfo.vertices(ii).input_filename_index = fscanf(fid,'%i',[1 1]);
	BRinfo.vertices(ii).type = fscanf(fid,'%i',[1 1]);
end
fclose(fid);


end




