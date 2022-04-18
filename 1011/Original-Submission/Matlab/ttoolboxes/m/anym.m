function A = anym( A , dir)
% B = anym(A , [dir])
% Recursive call of any.
% anym(c,[x1 x2 x3 ...  xn]) = any(...(any(any(any(c,x1),x2),x3),...),xn)
% If dir is not given, then any is called recursevily until the outcome is 0 or 1.
%
% Input:
%   c       the sequence
%   dir     the directions given as a vector. If the vector is empty, then the return value is (c~=0).
%
% Output 
%   c       any(...(any(any(any(c,x1),x2),x3),...),xn)
% 
% E.g.: anym([1 1 2; 0 1 1; 1 1 1])
%
% See also: any, allm
%
% Written by: tommsch, 2017

if(nargin==1)
    while(true);
        A=any(A);
        if(isscalar(A)); 
            break; end;
    end
elseif(nargin>1)
    if(isempty(dir)); 
        A=(A~=0); end;
    for i=dir
        A=any(A,i); end;
else
    A=false;
end

end