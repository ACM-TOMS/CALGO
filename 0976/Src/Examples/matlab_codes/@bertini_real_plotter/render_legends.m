
function br_plotter = render_legends(br_plotter)

f = fieldnames(br_plotter.legend);

a = [];
b = [];


for ii = 1:length(f)
	a = [a br_plotter.legend.(f{ii}).handles]; 
	b = [b br_plotter.legend.(f{ii}).text];
end



br_plotter.handles.legend = legend(a,b); 
set(br_plotter.handles.legend,'Location','NorthEast','Interpreter','none','FontSize',br_plotter.options.fontsizes.legend);


end


