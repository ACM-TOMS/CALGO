function v = parameters(M)
%PARAMETERS Return a structure of the Schwarz-Christoffel map parameters.

%   Copyright 1998 by Toby Driscoll.
%   $Id: parameters.m,v 2.1 1998/05/10 04:09:41 tad Exp $

v.diskmap = M.diskmap;
v.rectpolygon = M.rectpolygon;
v.rectaffine = M.rectaffine;
v.prevertex = vertex(M.rectpolygon);
