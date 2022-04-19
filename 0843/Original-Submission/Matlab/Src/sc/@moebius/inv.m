function M2 = inv(M1)
%INV    Invert a Moebius transformation.

%   Copyright (c) 1998 by Toby Driscoll.
%   $Id: inv.m,v 1.1 1998/07/01 20:14:05 tad Exp $

M2 = moebius;
M2.source = M1.image;
M2.image = M1.source;
M2.coeff = [1 -1 -1 1].*M1.coeff([1 3 2 4]);