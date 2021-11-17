function M1 = hplmap(M)
%HPLMAP Convert Schwarz-Christoffel disk map to a map from the half-plane.

%   Copyright 1998 by Toby Driscoll.
%   $Id: hplmap.m,v 2.1 1998/05/10 04:11:45 tad Exp $

p = polygon(M);
[z1,c1] = disk2hp(vertex(p),angle(p)-1,M.prevertex,M.constant);
M1 = hplmap(p,scmapopt(M),z1,c1);


