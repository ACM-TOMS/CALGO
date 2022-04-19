function JSR=blockjsr(varargin); 
% [ JSR ] = blockjsr( j1, j2, j3, ... ); 
% Returns the JSR of block diagonal matrices, given the JSR of the blocks.
% More precisely: Returns the bound for the JSR of a block diagonal matrix whose blocks have JSR j1, j2, ... . 
%
% Input: 
%   ji      either a 1x2 vector, or a scalar
%
% E.g.: blockjsr([1 0],2)
%
% Written by tommsch, 2017

MIN=max(cellfun(@min,varargin));
MAX=max(cellfun(@max,varargin));

%test input values
for i=1:nargin
    if(numel(varargin{i})>2);
        error('Entry %i is not an interval nor a scalar.',i);
    end
end

if(searchincellarray(MAX,varargin,1)); 
    JSR=MAX; 
else;
    JSR=[MIN MAX];
end;

end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing.   