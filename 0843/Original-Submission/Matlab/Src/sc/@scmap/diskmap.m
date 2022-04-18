function M = diskmap(M)
%DISKMAP Convert generic Schwarz-Christoffel map object to disk map.
%   DISKMAP(M) creates a diskmap object based on the polygon and
%   options contained in M.
%   
%   See the DISKMAP class documentation.

%   Copyright 1998 by Toby Driscoll.
%   $Id: diskmap.m,v 2.1 1998/05/10 04:23:04 tad Exp $

M = diskmap(M.polygon,M.options);
