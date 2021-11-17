function cr = crossrat(w,Q)
%CROSSRAT Crossratios of a triangulated polygon.
%   CROSSRAT(W,Q) returns the N-3 crossratios of the N polygon vertices
%   W defined by the triangulation as given by the quadrilateral graph
%   in Q.
%       
%   See also CRTRIANG, CRCDT, QLGRAPH.

%   Copyright 1998 by Toby Driscoll.
%   $Id: crossrat.m,v 2.1 1998/05/10 04:31:11 tad Exp $

wql = Q.qlvert;				% get size of wql right
wql(:) = w(Q.qlvert);
cr = (((wql(2,:)-wql(1,:)).*(wql(4,:)-wql(3,:)))./...
  ((wql(3,:)-wql(2,:)).*(wql(1,:)-wql(4,:)))).';
