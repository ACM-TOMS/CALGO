function [Au, idx ,idx2] = uniquecell(A,varargin)
% [Au, idx ,idx2] = uniquecell( A, [options])
% Same behaviour as unique, but for cell arrays consisting of numbers
% Converts cell arrays to strings and compares strings
%
% E.g.: uniquecell({1 2 1 [1 2] [1 3] [1 2]})
%       uniquecell({1 2 1 [1 2] [1 3] [1 2]; 10 20 10 [1 2] [1 3] [10 20]})
%
% See also: unique

    B = cellfun(@(x) num2str(x(:)'),A,'UniformOutput',false);
    if nargout > 2
        [~,idx,idx2] = unique(B, varargin{:});
        Au = A(idx);
    else
        [~,idx] = unique(B, varargin{:});
        Au = A(idx);
    end
end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing.   