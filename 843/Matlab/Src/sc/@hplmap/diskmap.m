function M1 = diskmap(M)
%DISKMAP Convert Schwarz-Christoffel half-plane map to a map from the disk.

%   Copyright 1998 by Toby Driscoll.
%   $Id: diskmap.m,v 2.1 1998/05/10 04:16:25 tad Exp $

p = polygon(M);
[z1,c1] = hp2disk(vertex(p),angle(p)-1,M.prevertex,M.constant);
M1 = diskmap(p,scmapopt(M),z1,c1);


