function M = plus(M,a)
%   Add a constant to the image of a diskmap (i.e., translate image).

%   Copyright (c) 1998 by Toby Driscoll.
%   $Id: plus.m,v 1.1 1998/06/29 23:02:02 tad Exp $

% May need to swap arguments
if isa(M,'double') & isa(a,'diskmap')
  tmp = M;
  M = a;
  a = tmp;
end

if length(a)==1 & isa(a,'double')
  M.center = M.center + a;
  M.scmap = M.scmap + a;
else
  error('Addition is not defined for these operands.')
end
