function BMat = genBasisMat( varStructure )
%GENBASISMAT   Generate basis matrices of moment matrix.
%   

computeMomentIdx = varStructure.computeMomentIdx;
sizeblk = varStructure.sizeblk;
numSupports = size(varStructure.supports,1);

blkIdx = varStructure.blkIdx;
rowIdx = varStructure.rowIdx;
colIdx = varStructure.colIdx;

momentIdxLow = computeMomentIdx(blkIdx, rowIdx, colIdx);
momentIdxUpp = computeMomentIdx(blkIdx, colIdx, rowIdx);
momentIdx = [momentIdxLow; momentIdxUpp];

varIdx = [varStructure.momentic; varStructure.momentic];

BMat = spones(sparse(momentIdx, varIdx, 1, sum(sizeblk.^2), numSupports));
end
