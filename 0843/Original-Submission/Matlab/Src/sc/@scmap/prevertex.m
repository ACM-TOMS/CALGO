function z = prevertex(M)
%PREVERTEX Extract a vector of the prevertices of an S-C map.

%   Copyright 1998 by Toby Driscoll.
%   $Id: prevertex.m,v 2.2 1998/05/13 23:29:29 tad Exp $

tmp = parameters(M);
if strmatch('prevertex',fieldnames(tmp))
  z = tmp.prevertex;
else
  msg = sprintf('Prevertices not defined for map of class %s\n',class(M));
  error(msg)
end
