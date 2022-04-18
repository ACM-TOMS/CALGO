function C = triu(C,k)
% D = triu(C,k)
% Same behaviour as triu for matrices
%
% E.g.:  triu({[2 3],[1 2 0];[2], [2;3];2 , 3})
%
% See also: triu, cell/diag
%
% Written by tommsch, 2018

if(nargin==1); 
    k=0; end;
for i=1:size(C,2)
    [C{max(1,i-k+1):end,i}]=deal([]);
end


end