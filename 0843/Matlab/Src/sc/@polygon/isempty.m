function t = isempty(p)
%   Returns true if there are no vertices.

%   Copyright 1998 by Toby Driscoll.
%   $Id: isempty.m,v 2.1 1998/05/10 03:52:39 tad Exp $

t = isempty(p.vertex);
