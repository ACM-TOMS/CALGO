%this function generates a filename, based on highest number already in
%directory.  

%could potentially be used elsewhere for autogeneration of
%names.

% daniel brake
% colorado state university, north carolina state university, notre dame
% mathematics and applied mathematics
% 2013-14
% danielthebrake@gmail.com



function newname = increment_name(basename)

file_list = dir([basename '_*']);
if ~isempty(file_list)
	filenumbers = zeros(1,length(file_list));
	for ii = 1:length(file_list)
		potential_number = str2double(file_list(ii).name(length(basename)+2:end-4));
		if isnan(potential_number)
			potential_number = -1;
		end
		
		filenumbers(ii) = potential_number;%subtract 4 to remove extension
	end
else
	filenumbers = 0;
end

newname = sprintf('%s_%i',basename,max(filenumbers)+1);


end
