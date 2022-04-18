function wp = eval(Md,zp)
%EVAL Evaluate differentiated SC map.
%   EVAL(MD,ZP), where MD is an SCMAPDIFF object and ZP is a vector of
%   points in canonical domain of the map, returns the derivative of the
%   map used to create MD at the points ZP.
%   
%   See also SCMAPDIFF, SCMAPDIFF/SUBSREF.

%   Copyright 1998 by Toby Driscoll.
%   $Id: eval.m,v 2.1 1998/05/10 04:25:00 tad Exp $

wp = evaldiff(Md.themap,zp);
