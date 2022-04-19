function ReadLOPInstance(filename)
% ReadLOPInstance(filename)
% ReadLOPInstance: Read LOP problem instance
% Input:
% filename: The file name which is to be read
%
% Global variable:
% LOPInstance:  Global variable LOPInstance, which have the read instance
%
% Created version 02/21/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 


global LOPInstance;

%Read the instance file
fid = fopen(filename);

%Read first line
tline = fgetl(fid);
[n,count] = sscanf(tline,'%d',1); % Matrix size
%If the file don't have the needed structure will give an error
if(count ~= 1)
	error('Wrong file format');
end

%Read the matrix values
% fscanf fills the array in column order
[A, count] = fscanf(fid,'%d',[n,n]);

%transpose the results
matrix = A';

%Store the values of the instance in a global variable
LOPInstance = {matrix, n};



