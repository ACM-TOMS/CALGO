function A = maxm(A, dir);
% B = maxm(A, [x1 x2 x3 ... xn])= max(...(max(max(all(A,x1),x2),x3),...),xn)
% If dir is not given, then all is called recursevily until the outcome is a scalar
%
% Input:
%   A       the array
%   dir     the directions given as a vector. If the vector is empty, then the return value is (c~=0).
%
% Output:
%   B       max(...(max(max(max(A,x1),x2),x3),...),xn)
%
% E.g.: maxm([1 2 3; 2 3 4])
%
% See also: minm, max
%
% Written by: tommsch, 2017

if(nargin==0 || isempty(A)); A=-inf; return; end;

if(nargin==1)
    while(true);
        A=max(A);
        if(isscalar(A)); 
            break; end;
    end
else
    if(isempty(dir)); 
        A=(A~=0); end;
    for i=dir
        A=max(A,[],i);
    end
end

end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing.   