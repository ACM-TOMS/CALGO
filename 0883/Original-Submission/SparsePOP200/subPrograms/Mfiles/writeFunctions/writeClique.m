function writeClique(fileId,clique)
fprintf(fileId,'# Maximal clique structure induced from csp graph\n');
% 2005 03 26
% H.Waki
% field 'NoC', 'maxC' and 'minC' of the structure 'clique' are sparse
% inputs. So, fprint can not output these. 
%
% I converted these into full input 
%

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
fprintf(fileId,'#Cliques=%d, maxC=%d, minC=%d\n',full(clique.NoC),full(clique.maxC),full(clique.minC));
[rowSize,colSize] = size(clique.Set);
for i=1:rowSize
	fprintf(fileId,'clique %2d : ',i);
	for j=1:colSize
		fprintf(fileId,'%3d',clique.Set(i,j));
	end
	fprintf(fileId,'\n'); 
end
fprintf(fileId,'\n'); 
return

% $Header: /home/waki9/CVS_DB/SparsePOPdev/subPrograms/writeFunctions/writeClique.m,v 1.1.1.1 2007/01/11 11:31:50 waki9 Exp $
