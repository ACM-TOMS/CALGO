function [found, ioo] = searchincellarray(newo,oo, vector_flag, len); 
% [ found, idx ] = searchincellarray(newo, oo, vector_flag, [len])
% Searches in cell array. More precisely:
% Searches for the first occurence of newo in the cell array oo.
%
% Input:
%   oo              cell array
%   newo            type of the elements in oo
%   vector_flag     If set to 1, then the elements in newo and oo are considered as a lineary collection of elements (vectors).
%                     This should be much faster.
%   len             Optional. length/size of the entries in oo
%                     If not given, these are computed by: len=cellfun(@(x) length(x(:)), oo); respectively: len=cellfun(@(x) size(x), oo,'UniformOutput',0);
%
% Output:
%   found           flag. 1 if someting is found, 0 if nothing found
%   idx             linear index of where it is found   
%
% E.g.: [found, ioo]=searchincellarray([1 2 3 4],{[1 3; 2 4],[1 2 3 4],[2 3 1]',[5]}, 0) %found=1, ioo=2
%       [found, ioo]=searchincellarray([1 2 3 4],{[1 3; 2 4],[1 2 3 4],[2 3 1]',[5]}, 1) %found=1, ioo=1
%
% Written by: tommsch, 2018

%#ok<*ALIGN>

found=false;
if(isempty(oo));
    ioo=0; 
    return; end;

if(vector_flag==1)
    newo=newo(:);
    lnewo=length(newo);
    if(nargin<=3); 
        len=cellfun(@(x) length(x(:)), oo); end;
    samelengthidx=find(len==lnewo);
    for i=1:numel(samelengthidx)
        val=isequal(newo,oo{samelengthidx(i)}(:));
        if(val);  
            ioo=samelengthidx(i);
            found=true;
            break; end;
    end
    
else
    szenewo=size(newo);
    if(nargin<=3); 
        len=cellfun(@(x) size(x), oo,'UniformOutput',0); end;
    samesizeidx=cellfun(@(x) isequal(szenewo,x),len);
    idx=cellfun(@(x) isequal(newo,x),oo(samesizeidx));
    ioo=find(samesizeidx); ioo=ioo(idx);
    found=any(idx);
end
if(found==0); 
    ioo=0; end;
end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing.   