function br_plotter = sphere_plot(br_plotter)


center = real(br_plotter.BRinfo.center(br_plotter.indices));
radius = real(br_plotter.BRinfo.radius);

switch length(br_plotter.indices)

	case 2
		num_samp = 100;
		to_2_pi = linspace(0,2*pi,num_samp);
		x = radius*cos(to_2_pi) + center(1);
		y = radius*sin(to_2_pi) + center(2);
		br_plotter.handles.sphere = plot(x,y,':k','LineWidth',2);
		
	case 3
		num_samp = 20;
		[x,y,z] = sphere(num_samp);
		h = patch(surf2patch(radius*x+center(1),radius*y+center(2),radius*z+center(3)));
		set(h,'FaceAlpha',0.1,'EdgeColor','none','FaceColor','k');
		br_plotter.handles.sphere = h;
	otherwise
		display('bad length(ind) in sphere_plot.')
		return;
end

switch br_plotter.dimension
	case 1
		% i have no idea right now how to tell from the data i already
		% collect whether this is bounded or not.
			br_plotter.is_bounded = 0;
		
	case 2
		if all( and(br_plotter.BRinfo.sphere_curve.edges(:,1)==br_plotter.BRinfo.sphere_curve.edges(:,2),br_plotter.BRinfo.sphere_curve.edges(:,3)==br_plotter.BRinfo.sphere_curve.edges(:,2)))
			br_plotter.is_bounded = 1;
		else
			br_plotter.is_bounded = 0;
		end
	otherwise
		
end



end
