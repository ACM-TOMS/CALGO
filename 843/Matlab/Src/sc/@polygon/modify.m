function [p,indx] = modify(p)
%MODIFY Modify a polygon graphically.
%   See MODPOLY for usage instructions.

%   Copyright 1998 by Toby Driscoll.
%   $Id: modify.m,v 2.1 1998/05/10 03:54:14 tad Exp $

[w,beta,indx] = modpoly(vertex(p),angle(p)-1);
p = polygon(w,beta+1);
