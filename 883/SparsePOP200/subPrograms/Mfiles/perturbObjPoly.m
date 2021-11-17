function [objPoly] = perturbObjPoly(objPoly,epsObj,randSeed)
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
rand('state',randSeed); 
nDim = objPoly.dimVar;  
mDim = objPoly.noTerms; 
objPoly.supports = [objPoly.supports; speye(nDim)];
objPoly.coef = [objPoly.coef; epsObj * (2*rand(nDim,1)-1)]; 
%epsObj*(2*rand(nDim,1)-1)
objPoly.noTerms = mDim + nDim; 
objPoly = simplifyPolynomial(objPoly);
return

% $Header: /home/waki9/CVS_DB/SparsePOPdev/subPrograms/perturbObjPoly.m,v 1.1.1.1 2007/01/11 11:31:50 waki9 Exp $
