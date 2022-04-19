function acc = accuracy(M)
%ACCURACY Apparent accuracy of Schwarz-Christoffel CR rectified map.

%   Copyright 1998 by Toby Driscoll.
%   $Id: accuracy.m,v 2.1 1998/05/10 04:08:27 tad Exp $

% Just return accuracy of underlying crdiskmap.
acc = accuracy(M.diskmap);
