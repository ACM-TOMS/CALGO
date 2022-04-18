function ReadPFSPInstance(filename)
% ReadPFSPInstance(filename)
% ReadPFSPInstance: Read PFSP problem instance
% Input:
% filename: The file name which is to be read
%
% Global variable:
% PFSPInstance:  Global variable PFSPInstance, which have the read instance
%
% Last version 12/18/2013. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 


global PFSPInstance;

filename
%Read the instance file
fid = fopen(filename);

%Read first line
tline = fgetl(fid);
%If the file don't have the needed structure will give an error
if(isempty(strfind(tline,'jobs')))
	error('Wrong file format');
end

%Read the variables
[A, count] = fscanf(fid,'%d',[1,5]);
v=num2cell(A);

%Store the readed values in variables
[jobs, machines, initial_seed, upper_bound, lower_bound] = v{:};

%Read end of line
tline = fgetl(fid);
%Read the matrix indicator line
tline = fgetl(fid);

if(isempty(strfind(tline,'processing times :')))
	error('Wrong file format');
end
%Read the matrix values
% fscanf fills the array in column order
[B, count] = fscanf(fid,'%d',[jobs,machines]);

%transpose the results
processingtimes = B';

%Store the values of the instance in a global variable
PFSPInstance = {processingtimes, machines, jobs};

%M = dlmread(filename, delimiter);



