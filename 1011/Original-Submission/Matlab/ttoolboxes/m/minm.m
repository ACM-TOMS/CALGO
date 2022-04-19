function A = minm(A, dir);
% B = minm(A, [x1 x2 x3 ... xn])= min(...(min(min(all(A,x1),x2),x3),...),xn)
% If dir is not given, then all is called recursevily until the outcome is a scalar
%
% Input:
%   A       the array
%   dir     the directions given as a vector. If the vector is empty, then the return value is (c~=0).
%
% Output:
%   B       min(...(min(min(min(A,x1),x2),x3),...),xn)
%
% E.g.: minm([1 2 3; 2 3 4])
%
% See also: maxm, min
%
% Written by: tommsch, 2017

 %#ok<*ALIGN>

if(nargin==0 || isempty(A)); 
    A=inf; 
    return; end;

if(nargin==1)
    while(true);
        A=min(A);
        if(isscalar(A)); 
            break; end;
    end
else 
    if(isempty(dir)); 
        A=(A~=0); end;
    for i=dir
        A=min(A,[],i);
    end
end

end