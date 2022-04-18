function M = mtimes(M,c)
%   Scale the image of a map by a complex constant.

%   Copyright (c) 1998 by Toby Driscoll.
%   $Id: mtimes.m,v 1.1 1998/06/29 23:09:57 tad Exp $

% May need to swap arguments
if isa(M,'double') & isa(c,'crdiskmap')
  tmp = M;
  M = c;
  c = tmp;
end

if length(c)==1 & isa(c,'double')
  M.affine = c*M.affine;
  M.scmap = c*M.scmap;
  M = center(M,c*M.center{1});
else
  error('Multiplication is not defined for these operands.')
end
