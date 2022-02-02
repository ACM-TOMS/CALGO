function writeSDPAformatData(fileId,param,SDPA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output SDPA sparse format file of SDP relaxation problem.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
lenBvect = size(SDPA.bVect,2);
[rowSize,colSize] = size(SDPA.sparseMatrix); 
fprintf('# Output : SDPA sparse format data\n'); 
fprintf('  File name = %s\n',param.sdpaDataFile);
fprintf('  mDim = %2d nBlock = %d\n',SDPA.mDim,SDPA.nBlock); 
fprintf('  size of bVect = 1 * %d\n',lenBvect);
fprintf('  size of sparseMatrix = %d * %d\n',rowSize,colSize);

fprintf(fileId, '* SDPA sparse format data\n');
fprintf(fileId, '* File name = %s\n',param.sdpaDataFile);
fprintf(fileId, '* mDim = %d, nBlock = %d\n',SDPA.mDim,SDPA.nBlock); 
fprintf(fileId, '* size of bVect = 1 * %2d\n',lenBvect);
fprintf(fileId, '* size of sparseMatrix = %2d * %2d\n',rowSize,colSize);
fprintf(fileId, '%d\n %d\n',SDPA.mDim,SDPA.nBlock); 
lenBolckStruct = length(SDPA.blockStruct); 
for j=1:lenBolckStruct
  fprintf(fileId,'%d ',SDPA.blockStruct(j));
end
fprintf(fileId, '\n');
  
for j=1:lenBvect
  fprintf(fileId,'%+15.10f ',SDPA.bVect(1,j));
end
fprintf(fileId,'\n'); 

for i=1:rowSize
  fprintf(fileId,'%6d ',SDPA.sparseMatrix(i,1));
  fprintf(fileId,'%6d ',SDPA.sparseMatrix(i,2)); 
  fprintf(fileId,'%6d ',SDPA.sparseMatrix(i,3)); 
  fprintf(fileId,'%6d ',SDPA.sparseMatrix(i,4));
  fprintf(fileId,'%+15.10f\n',SDPA.sparseMatrix(i,5));
end

return

% $Header: /home/waki9/CVS_DB/SparsePOPdev/subPrograms/writeFunctions/writeSDPAformatData.m,v 1.1.1.1 2007/01/11 11:31:51 waki9 Exp $
