function [ineqPolySys] = relax1EqTo2Ineqs(objPoly,ineqPolySys,eqTolerance)

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
if ~isempty(ineqPolySys)
	if iscell(ineqPolySys)
		nDim = ineqPolySys{1}.dimVar;  
	else
		nDim = ineqPolySys.dimVar;
	end
else
	nDim = objPoly.dimVar;
end
mDim = size(ineqPolySys,2);
pointer = mDim; 
for i=1:mDim
	if ineqPolySys{i}.typeCone == -1
		ineqPolySys{i}.typeCone = 1; 
		pointer = pointer + 1;
		ineqPolySys{pointer}.typeCone = 1; 
		ineqPolySys{pointer}.sizeCone = 1; 
		ineqPolySys{pointer}.degree = ineqPolySys{i}.degree; 
		ineqPolySys{pointer}.dimVar = ineqPolySys{i}.dimVar; 
		ineqPolySys{pointer}.noTerms = ineqPolySys{i}.noTerms; 
		ineqPolySys{pointer}.supports = ineqPolySys{i}.supports; 
		ineqPolySys{pointer}.coef = -ineqPolySys{i}.coef; 
		[I] = find(ineqPolySys{i}.supports * ones(nDim,1) == 0);
		if isempty(I)
			ineqPolySys{pointer}.noTerms = ineqPolySys{pointer}.noTerms + 1;
			ineqPolySys{pointer}.supports = [sparse(zeros(1,nDim)); ineqPolySys{pointer}.supports];
			ineqPolySys{pointer}.coef = [eqTolerance; ineqPolySys{pointer}.coef]; 
		else
			i = I(1); 
			ineqPolySys{pointer}.coef(i,1) =  ineqPolySys{pointer}.coef(i,1) + eqTolerance;
		end
	end
end
% size(ineqPolySys,2)

return 

% $Header: /home/waki9/CVS_DB/SparsePOPdev/subPrograms/relax1EqTo2Ineqs.m,v 1.1.1.1 2007/01/11 11:31:50 waki9 Exp $
