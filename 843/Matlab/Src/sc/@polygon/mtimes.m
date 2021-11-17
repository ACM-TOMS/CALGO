function r = mtimes(p,q)
%   Multiplication of a polygon by a scalar.

%   Copyright 1998 by Toby Driscoll.
%   $Id: mtimes.m,v 2.1 1998/05/10 03:54:47 tad Exp $

if isa(q,'polygon')
  if isa(p,'polygon')
    error('Function ''*'' not defined for two polygon objects.')
  end
  tmp = p;
  p = q;
  q = tmp;
end

r = p;
r.vertex = r.vertex*q;
