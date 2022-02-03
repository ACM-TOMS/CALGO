
function br_plotter = plot_projection(br_plotter)

ind = br_plotter.indices;

curr_axes = br_plotter.axes.main;
br_plotter.handles.projection = [];

fontsize = br_plotter.options.fontsizes.labels;


pi = br_plotter.BRinfo.pi;

labeltext = {};
br_plotter.handles.projection = zeros(br_plotter.BRinfo.dimension,1);
for ii = 1:br_plotter.BRinfo.dimension
	switch length(ind)
		case {2}
			k = plot([0 pi(ind(1),ii)],[0 pi(ind(2),ii)],'--k','Parent',curr_axes);
			
		case {3}
			k = plot3([0 pi(ind(1),ii)],[0 pi(ind(2),ii)],[0 pi(ind(3),ii)],'--k','Parent',curr_axes);
	end
	br_plotter.handles.projection(ii) = k;
	labeltext{ii} = ['\pi_' num2str(ii-1)];
end


switch length(ind)
	case {2}
		texthandles = text(pi(ind(1),:),pi(ind(2),:),labeltext,'HorizontalAlignment','right','FontSize',fontsize,...
			'Parent',curr_axes,'Color','r');
	case {3}
		texthandles = text(pi(ind(1),:),pi(ind(2),:),pi(ind(3),:),labeltext,'HorizontalAlignment','right','FontSize',fontsize,...
			'Parent',curr_axes,'Color','r');
end

br_plotter.handles.projection = [br_plotter.handles.projection; texthandles];
end



