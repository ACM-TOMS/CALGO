function A = allm( A ,dir)
% B = allm(A, [dir])
% Recursive call of all.
% allm(A,[x1 x2 x3 ...  xn]) = all(...(all(all(all(A,x1),x2),x3),...),xn)
% If dir is not given, then all is called recursevily until the outcome is 0 or 1.
%
% Input:
%   A       the array
%   dir     the directions given as a vector. If the vector is empty, then the return value is (c~=0).
%
% Output 
%   B       all(...(all(all(all(c,x1),x2),x3),...),xn)
% 
% E.g.: allm([1 1 2; 0 1 1; 1 1 1])
%
% See also: all, anym
%
% Written by: tommsch, 2017


if(nargin==1)
    while(true);
        A=all(A);
        if(isscalar(A)); 
            break; end;
    end
elseif(nargin>1) 
    if(isempty(dir)); 
        A=(A~=0); end;
    for i=dir
        A=all(A,i); end;
else
    A=true;
end

end