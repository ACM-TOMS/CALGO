function setupfig(br_plotter,varargin)
			
	if mod(length(varargin),2)~=0
		error('setupfig takes only pairs of arguments');
	end

% 	if ~isempty(which('getmondim'))
% 		mondim = getmondim;
% 		mondim(3:4) = min(mondim(3:4));
% 		br_plotter.window = [20 20 mondim(3)-60 mondim(4)-130];
% 	end

	fig = gcf;

	set(fig,'PaperPositionMode','auto');



	set(fig,'Name',sprintf('Bertini_real -- %s',br_plotter.filename),...
								'NumberTitle','off',...
								'Position',br_plotter.window);

	set(fig,'color','w');
	
	colormap(fig,br_plotter.options.colormap(br_plotter.options.num_colors))
	
	if isprop(fig,'SizeChangedFcn')
		set(fig,'SizeChangedFcn',@br_plotter.resizeui);
	else
		set(fig,'ResizeFcn',@br_plotter.resizeui);
	end
	br_plotter.figures.main = fig;
end %re: setupfig
