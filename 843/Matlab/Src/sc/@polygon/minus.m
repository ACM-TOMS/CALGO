function r = minus(p,q)
%   Translate a polygon by a scalar.

%   Copyright 1999 by Toby Driscoll.
%   $Id: minus.m,v 2.1 1999/09/03 23:25:27 tad Exp $

if isa(q,'polygon')
  if isa(p,'polygon')
    error('Function ''-'' not defined for two polygon objects.')
  end
  tmp = p;
  p = q;
  q = tmp;
end

if ~isa(q,'double') | length(q) > 1
  error('Only scalars can be added to a polygon.')
end
  
r = p;
r.vertex = r.vertex - q;
