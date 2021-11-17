function [vec] = ordering2vector(oo, l)
% [vec] = ordering2vector(oo, l)
% Takes an ordering and returns it as a vector of certain length.
%
% Input:
%   oo              ordering ( {[nxN], [nxM]} )
%   l               length of the output sequence, default: 100, but at least three repetitions of the periodic part
%
% Output: 
%   vec             ordering as a vector of length l
%
% Info:
%   If periodic part of oo is empty, and l is larger than size(oo,2), then a vector of size(oo,2) is returned only.
%   If the two matrices in oo have a different number of rows, the behaviour is undefined.
% 
% E.g.: ordering2vector({[2 2 3; 1 1 5],[4 4 4 0; 4 4 4 0]})
%       ordering2vector({ -1:-1:-40 , 1:50 })
%
% See also: ordering

l1=100;
l2=size(oo{1},2)+3*size(oo{2},2);
if(nargin==1); 
    l=max(l1,l2); end;

if(~isempty(oo{1}) && ~isempty(oo{2}));
    vec= [oo{1} repmat(oo{2},[1 ceil(l/size(oo{2},2))])];
elseif(isempty(oo{1}) && ~isempty(oo{2}));
    vec= repmat(oo{2},[1 ceil(l/size(oo{2},2))]);
elseif(~isempty(oo{1}) && isempty(oo{2}));
    vec=oo{1};
    loo=size(oo{1},2);
    if(loo<l); 
        l=loo; end;
else
    vec=zeros(size(oo{1},1),0);
    l=0;
end

vec=vec(:,1:l);

end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 