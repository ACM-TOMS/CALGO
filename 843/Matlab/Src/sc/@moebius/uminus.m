function M = uminus(M)
%UMINUS Negate a Moebius transformation.
%   Copyright (c) 1998 by Toby Driscoll.
%   $Id: uminus.m,v 1.1 1998/07/01 20:14:41 tad Exp $

M.coeff(1:2) = -M.coeff(1:2);
