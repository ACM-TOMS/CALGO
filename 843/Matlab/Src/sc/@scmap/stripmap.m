function M = stripmap(M)
%STRIPMAP Convert generic Schwarz-Christoffel map object to strip map.
%   STRIPMAP(M) creates a stripmap object based on the polygon and
%   options contained in M.
%   
%   See the STRIPMAP class documentation.

%   Copyright 1998 by Toby Driscoll.
%   $Id: stripmap.m,v 2.1 1998/05/10 04:24:34 tad Exp $

M = stripmap(M.polygon,M.options);
