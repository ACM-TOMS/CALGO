function D=supp(a, dim, amin)
% D = supp(a, dim, amin)
% Computes the support of a mask.
%
% Input: 
%   a       dim x N matrix
%   dim     the dimension of the array a (is needed since matlab deletes singleton dimensions)
%   amin    dim x 1 vector. Indices of the first entry in a
% 
% Output:
%   D       coordinates as column vectors
%   
% Note:
%   This function is much slower than the inverse function characteristic(). Use 'characteristic' if possible.
%
% E.g.: D=supp([1 1; 0 1],2,[1;-1])
%
% See also: characteristic
%
% Written by: tommsch, 2017


if(nargin~=3);
    error('Wrong arguments given.'); 
end;

a=(a~=0);
if(isa(a,'sym')); %test if <a> is symbolic
    a=isAlways(a,'Unknown','true'); 
end; 
SUM=sum(a);
while(length(SUM)>1); 
    SUM=sum(SUM); 
end; 

D=zeros(dim,SUM);
SIZE=size(a);
CO=cell(1,dim);

j=1;
for i=1:numel(a)
    [CO{:}]=ind2sub(SIZE,i);
    if(a(i)~=0); 
        COadd=[CO{:}]'-1;
        D(:,j)=setplus(COadd,amin); 
        j=j+1;
    end
end

end