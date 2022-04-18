function flag = issym(x)
% [ flag ] = issym( x )
% Tests if an object is symbolic.
% flag=logical(strcmp(class(x),'sym')); %#ok<STISA> 
%
% E.g.: issym([sym('23.2') 0])
%
% Written by: tommsch, 2018
    
    flag=logical(strcmp(class(x),'sym')); %#ok<STISA>
end