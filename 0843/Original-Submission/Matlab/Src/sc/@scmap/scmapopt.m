function opt = scmapopt(M,varargin)
%SCMAPOPT Options structure for a Schwarz--Christoffel map object.
%   Same as the regular SCMAPOPT, but the first argument is the options
%   structure contained in the map M.

%   Copyright 1998 by Toby Driscoll.
%   $Id: scmapopt.m,v 2.1 1998/05/10 04:24:27 tad Exp $

opt = scmapopt(M.options,varargin{:});
 