function s = summ(X, dim, absflag)
% s = summ(X, dim, [options]) 
% Recursive call of sum.
%
% Input:
%   X       the array to be summed up
%   dim     specifies the direction in which to sum
%               scalar      summ behaves like the function sum
%               vector      summ iteratively sums up in the directions specified by dim
%               empty       sum of all elements
%
% Options:
%   'abs'   sums up absolute values
%
% Note:
%   If X is empty, s equals 1.
%
% Output:
%   s       the sum
%
% E.g.: summ([2 3; -4 -5],[],'abs')
%
% See also: sum
%
% Written by: tommsch, 2017

if(nargin==1 || nargin>=2 && isempty(dim)); 
    allflag=true; 
else; 
    allflag=false; 
end;
if(nargin==3 && isequal(absflag,'abs')); 
    absflag=true; 
else; 
    absflag=false; 
end;

if(absflag); 
    X=abs(X); 
end;

if(allflag); 
    X=X(:); 
    s=sum(X); 
else;
    for i=1:length(dim);
        X=sum(X,dim(i));
    end
    s=X;
end

end