function varargout = get(map,varargin)
%GET    Get map parameters.
%   [VAL1,VAL2,...] = GET(F,'PROP1','PROP2',...) returns the values of the
%   map F corresponding to the requested properties. Valid properties
%   are: 
%   
%       polygon, options, prevertex, constant, center

% Copyright 1999 by Toby Driscoll.
% $Id: get.m,v 2.2 1999/09/30 23:04:45 tad Exp $

for j = 1:length(varargin)
  switch lower(varargin{j}(1:min(3,length(varargin{j}))))
   case 'pol'
    varargin{j} = map.scmap.polygon;
   case 'opt'
    varargin{j} = map.scmap.options;
   case 'pre'
    varargin{j} = map.prevertex;
   case 'con'
    varargin{j} = map.constant;
   case 'cen'
    varargin{j} = center(map);
   otherwise
    warning(sprintf('Property ''%s'' not recognized.\n',varargin{j}))
    varargin{j} = [];
  end
end
