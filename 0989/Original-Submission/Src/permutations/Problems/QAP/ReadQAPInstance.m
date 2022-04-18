function ReadQAPInstance(filename)
% ReadQAPInstance(filename)
% ReadQAPInstance: Read QAP problem instance
% Input:
% filename: The file name which is to be read
%
% Global variable:
% LOPInstance:  Global variable QAPInstance, which have the read instance
%
% Created version 03/03/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 


global QAPInstance;

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
distance = A';


%Read the matrix values
% fscanf fills the array in column order
[B, count] = fscanf(fid,'%d',[n,n]);

%transpose the results
flow = B';

%Store the values of the instance in a global variable
QAPInstance = {distance, flow, n};



