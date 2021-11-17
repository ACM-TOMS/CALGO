function get_indices(br_plotter)

	if br_plotter.BRinfo.num_variables<=4 % assumes the number of variables includes the homvar...
		vanilla_set_ind(br_plotter)
	else
		complicated_set_ind(br_plotter)
	end

end


function vanilla_set_ind(br_plotter)
	br_plotter.indices = 1:br_plotter.BRinfo.num_variables-1;
end


function complicated_set_ind(br_plotter)

	indices_of_nonconst_cols = get_nonconstind(br_plotter);

	if length(indices_of_nonconst_cols)>=4
		br_plotter.indices = get_user_indices(indices_of_nonconst_cols,br_plotter.BRinfo);
	else
		br_plotter.indices = indices_of_nonconst_cols;
	end

end




function indices_of_nonconst_cols = get_nonconstind(br_plotter)

	num_pts = min(1000,br_plotter.BRinfo.num_vertices);
	
	tmpdata = zeros(num_pts,br_plotter.BRinfo.num_variables-1);
	
	for ii = 1:num_pts
		tmpdata(ii,:) = br_plotter.BRinfo.vertices(ii).point(1:br_plotter.BRinfo.num_variables-1);
	end
	
	indices_of_nonconst_cols = find_constant_vars(tmpdata);
end