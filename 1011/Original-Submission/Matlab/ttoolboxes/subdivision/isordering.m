function x = isordering(o)
% x = isordering(S)
%  Tests whether input is an ordering.
% Format of ordering: 1x2 - cell array, 
%       first cell is an dim x N matrix
%       second cell is an dim x M matrix
%
% E.g.: isordering({[3],[1 2 3]})
%       isordering({zeros(1,0),[1 2 3; 0 0 0]})
%
% See also: constructordering
%
% Written by: tommsch, 2018

x=false;
if(~iscell(o)); 
    return; end;
if(size(o,2)~=2); 
    return; end;
%if(size(o{1},2)~=size(o{2},2)); 
%return; end;
if(~isempty(o{1}));
    %if(~isnumeric(o{1})); 
    %    return; end;
    if(~iswholenumber(o{1})); 
        return; end;
    %if(~isrow(o{1})); 
    %return; end;
end
if(~isempty(o{2}));
    %if(~isnumeric(o{2})); 
    %    return; end;
    if(~iswholenumber(o{2})); 
        return; end;
    %if(~isrow(o{2})); 
    %return; end;
end


x=true;

end