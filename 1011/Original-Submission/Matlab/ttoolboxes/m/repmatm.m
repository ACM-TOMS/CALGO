function C = repmatm(varargin)
% C = repmatm(varargin)
% Consistent behaviour of repmatm, regards to multi-dimensional applications.
%
% E.g.: repmatm(10,2)
%
% See also: repcellm, repmat
%
% Written by: tommsch, 2018

if(nargin==2);
    C=repmat(varargin{1},[varargin{2} 1]);
else
    C=repmat(varargin{1},varargin{2:end});
end

end