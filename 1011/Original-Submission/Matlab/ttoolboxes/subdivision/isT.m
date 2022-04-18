function x = isT(T, flag)
% x = isT(T, [options])
% Tries to find out if T is a cell array of transition matrices operators. Tests only formal stuff.
%
% Input:
%   T           the variable to be tested
%
% Options:
%   'noflat'    tests furthermore if T is a non-flat cell array of transition matrices
%   'flat'      tests furthermore if T is a flat cell array of transition matrices
%
% E.g.: isT(transitionmatrix(getS('1_rand')))
%
% See also: transitionmatrix, isS
%
% Written by: tommsch, 2017

%#ok<*ALIGN>

x=0;
if(~iscell(T)); 
    return; end

if(nargin==2 && isequal(flag,'flat'));  
    flat=1; 
else; 
    flat=0; end;
if(nargin==2 && isequal(flag,'noflat')); 
    noflat=1; 
else; 
    noflat=0; end;

%check whether we have flat or non-flat cell array
if(iscell(T{1})); 
    isflatflag=0; 
    T=[T{:}];
else
    isflatflag=1;     
end;
if( (noflat && isflatflag) || (flat && ~isflatflag) ); 
    return; end;

J=size(T,2);
SZE=size(T{1});
if(diff(SZE)); 
    return; end;
for j=1:J %test whether all T have the same dimension
    if(~isnumeric(T{j})); 
        return; end;
    if(~ismatrix(T{j})); 
        return; end;
    if(~isequal(size(T{j}),SZE)); 
        return; end;    
end

x=1;

end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 