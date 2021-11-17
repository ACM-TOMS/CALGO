function M = mtimes(M1,M2)
%MTIMES Multiply Moebius transformation by a scalar.

%   Copyright (c) 1998 by Toby Driscoll.
%   $Id: mtimes.m,v 1.1 1998/07/01 20:14:27 tad Exp $

% Make the first one moebius.
if isa(M1,'double')
  tmp = M1;
  M1 = M2;
  M2 = tmp;
end

C = M1.coeff;
if isa(M2,'double') & length(M2)==1
  C(1:2) = M2*C(1:2);
  M = moebius;
  M.coeff = C;
else
  error('Multiplication not defined for these operands.')
end
