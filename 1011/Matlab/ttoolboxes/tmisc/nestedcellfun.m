function [ ret ] = nestedcellfun(varargin)
% ret = nestedcellfun(varargin)
% Wrapper function calling cellfun for each cell in a nested array
% [ret] = nestedcellfun(func, C1, C2, ..., [options] ) applies the function func(x1,x2,...) to the contents of each cell of cell array C1,C2,... .
%
% Input:
%   func                    function handle
%   C1,C2,                  cell arrays, nested cell arrays or matrices. All of the Ci must have the same topology.
%
% Options:
%   'UniformOutput',false   Then ret has the same topology as C1, otherwise (default) it is a matrix. Thus the default behaviour is the same as for cellfun.
%
% Output:
%   ret                     Either a scalar or something with the same topology as C1
%                           If 'UniformOutput' is <false> and at least one cell array is given, the output is a cell array
%   
% Eg: nestedcellfun(@(x) ndimsm(x),{{[1 2 1]},[1 2; 2 3]},'UniformOutput',false)
%
% See also: flatten, unflatten
%
% Written by: tommsch, 2017

% Changelog: tommsch,   2019-03-04,     Now also handles multiple arguments which are not cell arrays

% if(~iscell(varargin{2})); 
%     varargin{2}=varargin(2); 
%     nocell=1; 
% else; 
%     nocell=0; 
% end;

[UniformOutput,varargin]=parsem('UniformOutput',varargin,true);
handle=varargin{1}; varargin(1)=[]; %the function handle

j=1;
C={};
Cf={};
cellflag=0;
while(~isempty(varargin))
    C{j}=varargin{1}; %#ok<AGROW>
    if(~iscell(C{j})); 
        C{j}=C(j);  %#ok<AGROW>
    else
        cellflag=1;
    end;
    Cf{j}=flatten(C{j}); %#ok<AGROW>
    j=j+1;
    varargin(1)=[];
end

if(isequal(UniformOutput,0) || strcmp(UniformOutput,'false'))
    UNFLAT=cellfun(handle,Cf{:},varargin{j:end},'UniformOutput',false);
    ret=unflatten(UNFLAT,C{1});
else
    ret=cellfun(handle,Cf{:},varargin{j:end});
end

if(~cellflag && ~UniformOutput); 
    ret=ret{1}; 
end;

end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 