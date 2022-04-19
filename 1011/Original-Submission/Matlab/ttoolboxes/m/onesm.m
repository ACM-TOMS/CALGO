function d = onesm(varargin)
% d = onesm( a1, ..., an) 
% d = onesm([a1 ... an]) 
% Consistent behaviour of ones, with regards to multi-dimensional applications.
%
% Example: x = onesm(3);
%
% See also eyem, zerosm
%
% Written by: tommsch, 2018

if(nargin==0); 
    d = [];
elseif(nargin>1);
    d=ones(varargin{:});
elseif(nargin==1 && isscalar(varargin{1}));
    d=ones(varargin{1},1);
elseif(nargin==1 && isempty(varargin{1}))
    d=ones(0,0);
else
    d=ones(varargin{:});
end;


end