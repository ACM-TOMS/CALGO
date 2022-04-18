function M = minus(M1,M2)
%MINUS  Subtract a scalar from a Moebius map.

%   Copyright (c) 1998 by Toby Driscoll.
%   $Id: minus.m,v 1.1 1998/07/01 20:14:10 tad Exp $

M = M1 + (-M2);
