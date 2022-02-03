function zr = rectangle(M)
%RECTANGLE Return the corners of the rectangle in the fundamental domain.

%   Copyright 1998 by Toby Driscoll.
%   $Id: rectangle.m,v 2.1 1998/05/10 04:22:19 tad Exp $

zr = prevertex(M);
zr = zr(corners(M));
