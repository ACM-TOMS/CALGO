function names = get_names(num_vars)

fid = fopen('names.out','r');

names = cell(num_vars,1);
for ii = 1:num_vars
	names{ii} = fgetl(fid);
end

fclose(fid);


end