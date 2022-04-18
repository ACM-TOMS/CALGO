function M = extermap(M)
%EXTERMAP Convert generic Schwarz-Christoffel map object to exterior map.
%   EXTERMAP(M) creates a extermap object based on the polygon and
%   options contained in M.
%   
%   See the EXTERMAP class documentation.

%   Copyright 1998 by Toby Driscoll.
%   $Id: extermap.m,v 2.1 1998/05/10 04:23:11 tad Exp $

M = extermap(M.polygon,M.options);
