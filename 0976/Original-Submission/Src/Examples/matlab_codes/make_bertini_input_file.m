function make_bertini_input_file(filename)

if nargin==0
	filename='input';
end
	
fid = fopen(filename,'w');

fprintf(fid,'CONFIG\n\ntracktype: 0;\n\nEND;\n\nINPUT\n\nvariable_group ;\n\nfunction ;\n\n\nEND;');

fclose(fid);


end
