function v = parameters(M)
%PARAMETERS Return a structure of the Schwarz-Christoffel map parameters.

%   Copyright 1998 by Toby Driscoll.
%   $Id: parameters.m,v 2.1 1998/05/10 04:18:46 tad Exp $

v.prevertex = M.prevertex;
v.constant = M.constant;