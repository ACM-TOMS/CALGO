function C = repcell(varargin)
% C = repcell(varargin)
% Returns the replicated argument as a cell array.
% Uses the same interface as repmat.
%
% E.g.: repcell(2 ,[1 2])
%
% See also: repmat

C=repmat(varargin(1),varargin{2:end});

end