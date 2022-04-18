function C = nondiag(C,k)
% C = nondiag(C,k)
% Extracts non-diagonal parts of 2-arrays and 2-cells.
%
% Input:
%   C       matrix or cell-matrix
%   k       which diagonal shall be removed.
%           k=0 is the main diagonal, k>0 is above, K<0 is below
%           (experimental) k=+-inf are the first or last diagonal
% 
% Output:
%   C       Same as C, but k-th diagonal is removed
%
% E.g.: nondiag({[1 2] [2 3] [1];[1] [40] [100; 10];[1] [2 2 2] [3;1]},1)
% 
% See also: triu
%
% Written by: tommsch, 2017

% YY Write nondiagm, diagm

if(nargin==1); 
    k=0; end;
if(k==inf); 
    k=size(C,1)-1; end;
if(k==-inf); 
    k=-size(C,1)+1; end;
for i=1:min(size(C,2)-k,size(C,1));
    if(i<=0||i+k<=0); 
        continue; end;
    if(iscell(C));
        C{i,i+k}=[];
    else
        C(i,i+k)=0;
    end
    
end

end