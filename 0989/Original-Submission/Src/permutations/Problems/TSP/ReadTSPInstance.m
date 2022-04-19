function ReadTSPInstance(filename)
% ReadTSPInstance(filename)
% ReadTSPInstance: Read TSP problem instance
% Input:
% filename: The file name which is to be read
%
% Global variable:
% TSPInstance:  Global variable TSPInstance, which have the read instance
%
% Created version 03/25/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 
%
% Last version 05/13/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es)

global TSPInstance;

%Read the instance file
fid = fopen(filename);

tline = fgetl(fid);
[n,count] = sscanf(tline,'%d',1); % Matrix size
%If the file don't have the needed structure will give an error
if(isempty(strfind(tline,'NAME')) && count ~= 1)
	error('Wrong file format');
end

if(not (isempty(strfind(tline,'NAME')))) %The instance has a structure

    %Read the variables
    [A, count] = textscan(fid,'%s',5,'delimiter','\n');
    v=A{1};

    for i=(1:length(v))
    splited(i,:) = strsplit(': ',char(v(i)));
    end

    %Create a key-value pairs map see containers.Map in matlab help
    keySet = splited(:,1);
    valueSet = splited(:,2);

    instParam = containers.Map(keySet,valueSet);

    %If the coord type is not implemented, will give an error
    if ( ~strcmp(instParam('EDGE_WEIGHT_TYPE'), 'EXPLICIT'))
        error(['Error the EDGE_WEIGHT_TYPE ', instParam('EDGE_WEIGHT_TYPE'), ' is not implemented']);
    end

    %Read Coordinate section
    tline = fgetl(fid);

    %Store the readed values in variables
    cities = instParam('DIMENSION') ;
    coord_type = instParam('EDGE_WEIGHT_TYPE');


    [n,count] = sscanf(cities,'%d',1); % Matrix size
    %If the file don't have the needed structure will give an error
    if(count ~= 1)
        error('Wrong file format');
    end

    %Read the matrix values
    % fscanf fills the array in column order
    [A, count] = fscanf(fid,'%d',n*(n+1)/2);

    uppertri = triu(ones(n),0); %create upper triangular matrix
    uppertri(uppertri==1) = A; %fills the upper triangular matrix with the values readed

    %Calculate simetrical Matrix
    distance = uppertri + uppertri';

elseif(count == 1) %The instance is only the distance matrix

    %Read the matrix values
    % fscanf fills the array in column order
    [A, count] = fscanf(fid,'%d',[n,n]);

    %transpose the results
    distance = A';
    
    
else
    error('Wrong file format');
end;

%Store the values of the instance in a global variable
TSPInstance = {distance, n};



