function r = plus(p,q)
%   Translate a polygon by a scalar.

%   Copyright 1998 by Toby Driscoll.
%   $Id: plus.m,v 2.1 1998/05/10 03:56:03 tad Exp $

if isa(q,'polygon')
  if isa(p,'polygon')
    error('Function ''+'' not defined for two polygon objects.')
  end
  tmp = p;
  p = q;
  q = tmp;
end

if ~isa(q,'double') | length(q) > 1
  error('Only scalars can be added to a polygon.')
end
  
r = p;
r.vertex = r.vertex + q;
