function out = tif(truth,yes,no)
% Ternary if operator
% out = truth & yes : no ;
%
% E.g.: tif(true, 1, 2)
    
if(truth); 
    out=yes; 
else; 
    out= no; 
end;
    
end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing.   