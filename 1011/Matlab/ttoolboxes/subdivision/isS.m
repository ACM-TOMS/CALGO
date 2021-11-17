function x = isS(S)
% x = isS(S)
% tries to find out if 'S' is a cell array of subdivision operators. Tests only formal stuff.
% to test if 'S' consists of meaningful subdiv operator, use getS(S) or getS(S,'bigcheck').
%
% E.g.: isS(getS(2))
%
% See also: isT, getS
%
% Written by: tommsch, 2016

x=false;

if(~iscell(S)); 
    return; end;

sizeS=size(S,1); 

if(size(S,2)~=4); %tests if there are at least 4 elements per row
    return; end; 

for i=1:sizeS
    a=S{i,1};
    M=S{i,2};
    D=S{i,3};
    name=S{i,end};
    
    if(~issquare(M));  %test if dilation matrices are square
        return; end; 
    if(~iswholenumber(M)); %test if dilation matrices have integer values
        return; end; 
    
    if(~isa(a,'sequence')); %tests if masks are a cell array
        return; end; 
    
    if(~isnumeric(D)); %test if M is symbolic  
        return; end;       
    if(size(D,1) ~= size(M,1)); %tests if digits have the right dimension
        return; end; 
    
    if(~ischar(name));%test if name is a string
        return; end; 
end

x=true;

end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 