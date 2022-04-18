







function [indices] = get_user_indices(suggested_indices, BRinfo)

Title = 'Bertini_real variable plotting options';

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

names = make_names(suggested_indices, BRinfo);





Prompt(1,:) = {' x   ','var1',[]};
Formats(1,1).type = 'list';
Formats(1,1).style = 'popupmenu';
Formats(1,1).size = [80 0];
Formats(1,1).items = names;%{'Black','White','Red','Blue','Green','Yellow','Orange'};

default_answers(1).var1 = suggested_indices(1);



Prompt(2,:) = {' y   ','var2',[]};
Formats(2,1).type = 'list';
Formats(2,1).style = 'popupmenu';
Formats(2,1).size = [80 0];
Formats(2,1).items = names;%{'Black','White','Red','Blue','Green','Yellow','Orange'};

default_answers.var2 = suggested_indices(2);


if BRinfo.num_variables-1>=3
	
	names{end+1} = '-';
	Prompt(3,:) = {' z   ','var3',[]};
	Formats(3,1).type = 'list';
	Formats(3,1).style = 'popupmenu';
	Formats(3,1).size = [80 0];
	Formats(3,1).items = names;%{'Black','White','Red','Blue','Green','Yellow','Orange'};

	default_answers.var3 = suggested_indices(3);

end





[answer,user_cancelled] = inputsdlg(Prompt,Title,Formats,default_answers,Options);


if user_cancelled
	error('cancelled variable choice.  cannot continue');
end



indices(1) = answer.var1;
indices(2) = answer.var2;

if answer.var3 ~= BRinfo.num_variables
	indices(3) = answer.var3;
end


end





function names = make_names(suggested_indices, BRinfo)
	names = BRinfo.var_names;
	
	
	max_name_length = 0;
	
	for ii = 1:BRinfo.num_variables-1
		if length(names{ii}) > max_name_length
			max_name_length = length(names{ii});
		end
	end
	

	for ii = 1:BRinfo.num_variables-1
		base_length = length(names{ii});
		for jj = base_length:max_name_length
			names{ii} = sprintf('%s ',names{ii});
		end
		
		if isempty(find(suggested_indices==ii, 1))
			names{ii} = sprintf('%s (constant)',names{ii});
		end
	end


end











