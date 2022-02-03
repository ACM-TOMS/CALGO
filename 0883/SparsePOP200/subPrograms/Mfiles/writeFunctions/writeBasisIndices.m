function writeBasisIndices(fileId,basisIndices)
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
mDim = size(basisIndices,2); 
fprintf(fileId,'# basisIndices --- the set of nonzero coordinates\n'); 
fprintf(fileId,'  used for each equality and inequality\n'); 
for i=1:mDim
	cliqueSize = length(basisIndices{i});
	fprintf(fileId,'%3d:',i);
	for j=1:cliqueSize
		fprintf(fileId,'%3d',basisIndices{i}(j));
	end
	fprintf(fileId,'\n');
end
fprintf(fileId,'\n'); 
return

% $Header: /home/waki9/CVS_DB/SparsePOPdev/subPrograms/writeFunctions/writeBasisIndices.m,v 1.1.1.1 2007/01/11 11:31:50 waki9 Exp $
