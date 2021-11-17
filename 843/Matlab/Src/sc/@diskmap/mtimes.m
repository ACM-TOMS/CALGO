function M = mtimes(M,c)
%   Scale the image of a map by a complex constant.

%   Copyright (c) 1998 by Toby Driscoll.
%   $Id: mtimes.m,v 1.2 1998/06/29 23:01:00 tad Exp $

% May need to swap arguments
if isa(M,'double') & isa(c,'diskmap')
  tmp = M;
  M = c;
  c = tmp;
end

M.constant = c*M.constant;
M.scmap = c*M.scmap;
M.center = c*M.center;