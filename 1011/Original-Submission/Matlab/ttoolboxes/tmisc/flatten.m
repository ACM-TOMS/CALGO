function [ Cflat ] = flatten(C)
% Cflat = flatten(C)
% Changes a nested cell array to a flat cell array
% This function is slow, since it treats every cell entry by its own! This function should be rewritten.
%
% Input: 
%       C           a nested cell array of arbitrary topology
%                   If C is not a cell array, it will be a cell array afterwards containing of the one element C
% Output:
%       Cflat       the flat cell array
%
% Notes:
% This function is used in 'nestedcellfun' to apply a function to each element of a nested cell array.
% Or, when used together with unflatten(Cflat, C) something like a "foreach" loop can be made manually
%       Cflat=flatten(C)
%       for i=1:size(Cflat,2)
%           %do something with c{i}
%       end
%       C=unflatten(Cflat, C)
%
% E.g.: flatten({{2 3 {4; 5} 6} 7})
% 
% See also: unflatten, nestedcellfun
%
% Written by: tommsch, 2017

if(~iscell(C));  %#ok<ALIGN>
    Cflat{1}=C; 
    return; end;
Cflat={};
for i=1:numel(C)
    if(iscell(C{i}));
        Cflat=[Cflat flatten(C{i})]; %#ok<AGROW>
    else
        Cflat{end+1}=C{i}; %#ok<AGROW>
    end
end


end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 
