function array = savetocellarray(value, idx, array)
% array = savetocellarray(value, idx, array)
% Stores values in a cell array corresponding to a linear index-vector.
% More precisely: Stores <val> in the cell array <array>, corresponding to the places where idx=1.
%
%Input:
%   value       values to save with length nnz(idx), or a scalar
%   idx         logical array with length length([array{:}])
%   array       the cell array
%
% Output:
%   array       The array with the stored values
%
% E.g.: vdisp(savetocellarray([1 2 3],logical([1 1 0 0 1 ]),{[-1 -2], [-3 -4 -5]}))
%
% Written by: tommsch, 2018

    if(nnz(idx)==0); 
        return; end;
    if(isempty(idx)); 
        idx=true(size(array)); end;
    if(isscalar(value)); 
        value=repmat(value,1,nnz(idx)); end;
    L=cellfun(@length,array);
    
    idx=mat2cell(idx,1,L);
    value=mat2cell(value,1,cellfun(@nnz,idx));
    for i=1:size(array,2)
        array{i}(idx{i})=value{i};
    end


end

