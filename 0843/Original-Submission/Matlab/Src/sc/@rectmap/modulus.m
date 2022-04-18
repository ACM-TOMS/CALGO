function mu = modulus(M)
%MODULUS Conformal modulus of the generalized quadrilateral.
%   Returns the conformal modulus of the polygon in the rectmap (the
%   aspect ratio of the source rectangle).

%   Copyright 1998 by Toby Driscoll.
%   $Id: modulus.m,v 2.1 1998/05/10 04:21:49 tad Exp $

z = M.prevertex;
mu = max(imag(z)) / (2*max(real(z)));
