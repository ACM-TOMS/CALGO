function C = repcellm(varargin)
% C = repcellm(varargin)
% Consistent behaviour of repcellm, regards to multi-dimensional applications.
% Uses the same interface as repmatm.
%
% E.g.: repcellm(10 ,2)
%
% See also: repcell, repmatm, repmat
%
% Written by: tommsch, 2018

C=repmatm(varargin(1),varargin{2:end});

end