function [objPoly0,ineqPolySys0,lbd0,ubd0] = saveOriginalPOP(objPoly,ineqPolySys,lbd,ubd)

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is a component of SparsePOP 
% Copyright (C) 2007 SparsePOP Project
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if length(lbd) ~= length(ubd)
   error('length of lower bound ~= length of upper bound');
end
if size(lbd,1) ~=1
   error('Please define lower bound as row vector'); 
end
if size(ubd,1) ~=1
   error('Please define upper bound as row vector'); 
end

eps = 1.0e-10;
DiffIdx = find(ubd - lbd <= -eps);
if ~isempty(DiffIdx)
   error('lower and upper bound is inconsistent because lower bound is larger than upper at some index'); 
end

objPoly0 = objPoly;
ineqPolySys0 = ineqPolySys;
lbd0 = lbd;
ubd0 = ubd;

return

% $Header: /home/waki9/CVS_DB/SparsePOPdev/subPrograms/saveOriginalPOP.m,v 1.3 2007/01/12 09:20:43 waki9 Exp $
