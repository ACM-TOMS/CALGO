function d = zerosm(varargin)
% d = zerosm( a1, ..., an) 
% d = zerosm([a1 ... an]) 
% Consistent behaviour of zero, regards to multi-dimensional applications.
%
% Example: x = zerosm(3);
%
% See also eyem, onesm
%
% Written by: tommsch, 2018

if(nargin==0); 
    d = [];
elseif(nargin>1);
    d=zeros(varargin{:});
elseif(nargin==1 && isscalar(varargin{1}));
    d=zeros(varargin{1},1);
elseif(nargin==1 && isempty(varargin{1}))
    d=zeros(0,0);
else
    d=zeros(varargin{:});
end;


end