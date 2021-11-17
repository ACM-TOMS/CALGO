function Md = diff(M)
%DIFF   Differentiated SC map object.
%   DIFF(M) returns an object formally representing the derivative of
%   the map M.

%   Copyright 1998 by Toby Driscoll.
%   $Id: diff.m,v 2.1 1998/05/10 04:22:38 tad Exp $

Md = scmapdiff(M);
