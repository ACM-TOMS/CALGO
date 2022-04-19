function br_plotter = change_text_size(br_plotter,source,event)
			



[text_options, user_cancelled] = getopt(br_plotter);
			

if user_cancelled
	return;
end





set_legend_text_size(br_plotter,[],[],text_options.legend);


set_axis_text_size(br_plotter,[],[],text_options.axis);


set_label_text_size(br_plotter,[],[],text_options.labels);



			
end









%http://www.mathworks.com/matlabcentral/fileexchange/25862-inputsdlg--enhanced-input-dialog-box--v2-0-5-


function [answer, cancelled] = getopt(br_plotter)

Title = 'Bertini_real FontSize options';

%%%% SETTING DIALOG OPTIONS
% Options.WindowStyle = 'modal';
Options.Resize = 'on';
Options.Interpreter = 'none';
Options.CancelButton = 'on';
Options.ApplyButton = 'off';
Options.ButtonNames = {'Continue','Cancel'}; %<- default names, included here just for illustration


Prompt = {};
Formats = {};
default_answers = struct([]);





Prompt(1,:) = {'Choose your FontSize options',[],[]};
Formats(1,1).type = 'text';
Formats(1,1).size = [-1 0];
% Formats(1,1).span = [1 2]; % item is 1 field x 4 fields





Prompt(2,:) = {'Vertex Labels       ', 'labels',[]};
Formats(2,1).type = 'edit';
Formats(2,1).format = 'integer';
Formats(2,1).limits = [0 96]; % 9-digits (positive #)
Formats(2,1).size = 60;
default_answers(1).labels = br_plotter.options.fontsizes.labels;


Prompt(3,:) = {'Axis                       ', 'axis',[]};
Formats(3,1).type = 'edit';
Formats(3,1).format = 'integer';
Formats(3,1).limits = [0 96]; % 9-digits (positive #)
Formats(3,1).size = 60;
default_answers(1).axis = br_plotter.options.fontsizes.axis;

Prompt(4,:) = {'Legend                 ', 'legend',[]};
Formats(4,1).type = 'edit';
Formats(4,1).format = 'integer';
Formats(4,1).limits = [0 96]; % 9-digits (positive #)
Formats(4,1).size = 60;
default_answers(1).legend = br_plotter.options.fontsizes.legend;




[answer,cancelled] = inputsdlg(Prompt,Title,Formats,default_answers,Options);


end

