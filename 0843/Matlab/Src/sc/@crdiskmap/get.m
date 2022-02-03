function varargout = get(map,varargin)
%GET    Get map parameters.
%   [VAL1,VAL2,...] = GET(F,'PROP1','PROP2',...) returns the values of the
%   map F corresponding to the requested properties. Valid properties
%   are: 
%   
%       polygon, options, prevertex, crossratio, affine, qlgraph,
%       original, qdata, center

% Copyright 1999 by Toby Driscoll.
% $Id: get.m,v 1.3 1999/09/30 23:00:43 tad Exp $

for j = 1:length(varargin)
  switch lower(varargin{j}(1:min(3,length(varargin{j}))))
   case 'pol'
    varargin{j} = map.scmap.polygon;
   case 'opt'
    varargin{j} = map.scmap.options;
   case 'pre'
    param = parameters(map);
    varargin{j} = param.prevertex;
   case 'cro'
    varargin{j} = map.crossratio;
   case 'aff'
    varargin{j} = map.affine;
   case 'qlg'
    varargin{j} = map.qlgraph;
   case 'ori'
    varargin{j} = map.original;
   case 'qda'
    varargin{j} = map.qdata;
   case 'cen'
    varargin{j} = map.center{1};
   otherwise
    warning(sprintf('Property ''%s'' not recognized.\n',varargin{j}))
    varargin{j} = [];
  end
end
