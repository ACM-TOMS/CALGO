function momentMat = genLowTriMomentMat( idxSupports, I01 )
%GENLOWTRIMOMENTMAT   Generate lower triangular part of moment matrix
%   
%   INPUT: 
%      idxSupports:  cell 
%      I01        :  a logical vector which specifies 0-1 binary variables
%
%   OUTPUT:
%      momentMat  :  struct which contains the following fields
%         .blkIdx   : vector of indices of blocks (i.e., value of k).
%         .rowIdx   : vector of indices of row. (i)
%         .colIdx   : vector of indices of column. (j)
%         .supports : matrix whose l-th row represents the support of
%                the (rowIdx(l), colIdx(l)) element 
%                of the blkIdx(l)-th block of the moment matrix.
%         .computeMomentIdx: function handle to compute l from (k,i,j)

numblk = length(idxSupports);

sizeblk = zeros(1,numblk);
for ii = 1:numblk
    sizeblk(ii) = size(idxSupports{ii},1);
end
pointerBlk = cumsum([0,sizeblk]);
numelTrilBlk = (sizeblk+1).*sizeblk / 2;

eidx = cumsum(numelTrilBlk);
sidx = [0, eidx(1:end-1)] + 1;

momentMat.blkIdx = zeros(eidx(end), 1);
momentMat.rowIdx = zeros(eidx(end), 1);
momentMat.colIdx = zeros(eidx(end), 1);

rowIdx = cell(numblk,1);
colIdx = cell(numblk,1);
for kk = 1:numblk
    [mrow,mcol] = find(tril(true(sizeblk(kk)))); %Todo: use enumCombWithRep
    momentMat.blkIdx(sidx(kk):eidx(kk)) = kk;
    momentMat.rowIdx(sidx(kk):eidx(kk)) = mrow; 
    momentMat.colIdx(sidx(kk):eidx(kk)) = mcol; 
    rowIdx{kk} = pointerBlk(kk) + mrow;
    colIdx{kk} = pointerBlk(kk) + mcol;
end
rowIdx = cell2mat(rowIdx);
colIdx = cell2mat(colIdx);

tmpSupports = cell2mat(idxSupports);
momentMat.supports = tmpSupports(rowIdx,:) + tmpSupports(colIdx,:);
momentMat.supports(:,I01) = spones(momentMat.supports(:,I01));

sizeblk = sizeblk';
numelBlk = sizeblk.^2;
cumsumNumelBlk = cumsum([0;numelBlk]);
momentMat.computeMomentIdx = @(blk,row,col)...
    cumsumNumelBlk(blk) + (row-1) .* sizeblk(blk) + col;

end