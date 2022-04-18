function n = size(p,m)
%   Number of vertices.

%   Copyright 1998 by Toby Driscoll.
%   $Id: size.m,v 2.1 1998/05/10 04:00:22 tad Exp $

if nargin == 1
  n = [length(p.vertex) 1];
elseif m == 1
  n = length(p.vertex);
else
  n = 1;
end

  