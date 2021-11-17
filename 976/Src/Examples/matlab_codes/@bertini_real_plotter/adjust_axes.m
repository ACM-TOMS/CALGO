function br_plotter = adjust_axes(br_plotter)

curr_axes = br_plotter.axes.main;


lower = min(br_plotter.fv.vertices,[],1);
upper = max(br_plotter.fv.vertices,[],1);


ten_percent = 0.10 *(upper-lower);

ten_percent(ten_percent<0.01) = 0.01;


lower = lower - ten_percent;
upper = upper + ten_percent;

num_plotted_vars = size(br_plotter.fv.vertices,2);

if num_plotted_vars==3
    axis([lower(1) upper(1) lower(2) upper(2) lower(3) upper(3)]);
	
else
    axis([lower(1) upper(1) lower(2) upper(2)]);
end

widths = upper-lower;
ratios = [];
for ii = 1:num_plotted_vars
	for jj = ii+1:num_plotted_vars
		ratios(end+1) = min(widths(ii),widths(jj)) / max(widths(ii),widths(jj));
	end
end
min_ratio = max(ratios);
max_ratio = max(ratios);
if and( min_ratio > 0.5 , max_ratio < 2)
	set(curr_axes,'dataaspectratio',[1 1 1]);
end
% 

end
