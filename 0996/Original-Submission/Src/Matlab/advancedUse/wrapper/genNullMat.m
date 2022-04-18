function NMat = genNullMat( varStructure )
%GENNULLMAT Generate a matrix for the nullspace representation of a moment matrix

dimBlk = varStructure.sizeblk';
numelBlk = dimBlk.^2;
cumsumNumelBlk = cumsum([0;numelBlk]);
numSupports = size(varStructure.supports,1); % length(supSet01.ia);

function idx = computeMomentIdx(blk,row,col)
    idx = cumsumNumelBlk(blk) + (row-1) .* dimBlk(blk) + col;
end

blkIdx = varStructure.blkIdx;
rowIdx = varStructure.rowIdx;
colIdx = varStructure.colIdx;

momentIdxLow = computeMomentIdx(blkIdx, rowIdx, colIdx);
momentIdxUpp = computeMomentIdx(blkIdx, colIdx, rowIdx);

numEq = length(varStructure.momentic) - length(varStructure.ia);
workspace = repmat(1:numEq,4,1);
NMatrow = workspace(:);
workspace = [ones(2,numEq);-ones(2,numEq)];
val = workspace(:);

momentIdx = zeros(numEq*4,1);
pointer = 0;
for kk = 1:numSupports
    wherekk = varStructure.momentic == kk;
    workvec1 = momentIdxLow(wherekk)';
    workvec2 = momentIdxUpp(wherekk)';
    workspace = [workvec1(1:end-1);workvec2(1:end-1);...
        workvec1(2:end);workvec2(2:end)];
    numElForkk = 4*(sum(wherekk)-1);
    momentIdx( pointer + (1:numElForkk) ) = workspace(:);
    pointer = pointer + numElForkk;
end

NMat = sparse(NMatrow, momentIdx, val, numEq, sum(numelBlk));

end