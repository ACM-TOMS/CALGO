function [objPoly, I01, Icomp, relaxOrder, params] = ...
        genPOParrow(degree,a,b,c,l,isBin,addComplement)
%% [objPoly, I01, Icomp, relaxOrder, params]=genPOParrow(degree,a,b,c,l,isBin,addComplement)    
%% Input
%%   degree  : the degree of the polynomial objective function
%%   a       : the block size
%%   b       : the number of conseccutive variables 
%%   c       : the number of common variables over all the blocks
%%   l       : blockSize, the number of blocks
%%   isBin   : 1 for the binary constraint, 0 for the box constraint
%%   addComplement   : 1 for adding random complementarity constraints 
%%
%% a, b, c, l determines the Hessian matrix of the objective function;
%% the size of the Hessian matrix = the number of variabls is given by 
%%     n = (a - b) * l + b + c;
%% for example, if a = 3, l = 2, b = 1 and c = 2, then the Hessian matrix 
%% if 7 x 7 matorix of the form 
%%  ***00**
%%  ***00**
%%  *******
%%  00*****
%%  00*****
%%  *******
%%  *******
%% a, b, c and l must satisfy 
%% 0 <= b < a, 0 <= c, 2 <= l
%%
%% block-diagonal-type
%% if a >=1, b = 0, c = 0 and l >= 2. 
%% simple arrow-type
%% if a = 1, b = 0, c >= 1 and l >= 2. 
%% block-bandwidth-type
%% if 2 <= a, 1 <= b < a, c = 0 and l >= 2. 

    relaxOrder = ceil(degree/2); 
    params.delta1 = 1.0e-4;
    params.delta = 0;

    %%%%%%%%%%
    randSeed = 2017;
    n = (a - b) * l + b + c;
    I01 = isBin .* true(1, n);
    probData.nDim = n;
    probData.degree = degree;
    probData.randSeed = randSeed;
    probData.blockSize = a;
    probData.consecVarSize = b;
    probData.commonVarSize = c;
    probData.noBlocks = l;
    %%%%%%%%%%
    
    rng('default');
    rng(probData.randSeed);

    nDimBk = probData.blockSize ...
        + (probData.noBlocks - 1) * (probData.blockSize - probData.consecVarSize);
    objPoly.typeCone = 1;
    objPoly.sizeCone = 1;
    objPoly.dimVar = probData.nDim;
    objPoly.degree = probData.degree;
    probData.blkStartIdx=...
        [0:probData.noBlocks-1] * (probData.blockSize - probData.consecVarSize);
    probData.commonBlkIndices = (nDimBk+1):(nDimBk+probData.commonVarSize);
    binaryIdxAll = find(I01); %find(objPoly.typeVar == 1);
    nonBinaryIdxAll = find(~I01); %find(objPoly.typeVar == 0);
    supportsT = [];
    clique.Set = sparse(probData.noBlocks, probData.nDim);
    for p=1:probData.noBlocks
        oneBlock=[(probData.blkStartIdx(p)+1):(probData.blkStartIdx(p)...
            + probData.blockSize), probData.commonBlkIndices];
        clique.Set(p, :) = sparse(1, oneBlock, true, 1, probData.nDim);
        binaryIdx = intersect(oneBlock, binaryIdxAll);
        nonBinaryIdx = intersect(oneBlock, nonBinaryIdxAll);
        [oneSupportsT] = fullSup2(probData.nDim, binaryIdx, nonBinaryIdx, objPoly.degree);
        supportsT = [supportsT, oneSupportsT];
    end
    objPoly.supports = supportsT';
    objPoly.supports = unique(objPoly.supports, 'rows', 'first', 'legacy');
    objPoly.noTerms = size(objPoly.supports, 1);
    objPoly.coef = 2*rand(objPoly.noTerms, 1) - ones(objPoly.noTerms, 1);
    x0 = zeros(size(objPoly.supports, 2), 1);
    objPoly.UbdObjVal = evalPoly(objPoly, x0);

    %%%%%%%%%%
    Icomp = [];
    if addComplement
        Icomp = complementarityConstraint(2, 2, clique, randSeed);
    end
    %%%%%%%%%
end


function [supIdxVect] = fullSup2(nDim,binaryIdx,nonBinaryIdx,relaxOrder,supIdxVect1)
    if relaxOrder == 0
        supIdxVect = sparse(nDim,1); 
     elseif (relaxOrder == 1) 
         if (nargin == 4) || isempty(supIdxVect1)
             allIdx = union(nonBinaryIdx,binaryIdx);
             sDim = length(allIdx);
             supIdxVect = sparse(allIdx,(2:(sDim+1)),1,nDim,sDim+1,sDim+1);
         else
             supIdxVect = supIdxVect1;
         end
    elseif isempty(nonBinaryIdx)
        relaxOrder = min(relaxOrder,length(binaryIdx)); 
        I = [];
        J = [];
        colIdx = 1;
        for i=1:relaxOrder
            oneI = nchoosek(binaryIdx,i)'; 
            [kDim,mDim] = size(oneI); 
            oneI = reshape(oneI,kDim*mDim,1);
            oneJ = reshape(repmat(1:mDim,kDim,1),kDim*mDim,1);
            I = [I;oneI];
            J = [J;oneJ+colIdx];
            colIdx = colIdx + max(oneJ);
        end
        supIdxVect = sparse(I,J,1,nDim,colIdx,length(I));        
    else % (relaxOrder >= 1) && ~isempty(nonBinaryIdx)
        allIdx = union(nonBinaryIdx,binaryIdx);
        sDim = length(allIdx);
        supIdxVect1 = sparse(allIdx,(2:(sDim+1)), ... 
            1,nDim,sDim+1,sDim+1);
        mDim1 = size(supIdxVect1,2); 
        rOrderF = floor(relaxOrder/2);
        rOrderC =  ceil(relaxOrder/2); 
        if rOrderF >= 2
            supIdxVectF = ... 
                fullSup2(nDim,binaryIdx,nonBinaryIdx,... 
                rOrderF,supIdxVect1); 
        else
            supIdxVectF = supIdxVect1; 
        end
        mDimF = size(supIdxVectF,2);
        if rOrderC == rOrderF
            supIdxVectC = supIdxVectF;
            mDimC = mDimF; 
        else  % rOrderF < rOrderC          
            % supIdxSetF \times supIdxSet1 to create supIdxVectC
            mDim = mDim1*mDimF;         
            supIdxVectC = ... 
                reshape(repmat(supIdxVectF,mDim1,1),nDim,mDim);
            supIdxVect1 = ... 
                reshape(repmat(supIdxVect1,1,mDimF),nDim,mDim); 
            supIdxVectC = supIdxVectC + supIdxVect1;
            supIdxVectC(binaryIdx,:) = ... 
                spones(supIdxVectC(binaryIdx,:)); 
            supIdxVectC = ... 
                unique(supIdxVectC','rows','first','legacy')'; 
            mDimC = size(supIdxVectC,2);           
        end
        mDim = mDimF*mDimC; 
        % supIdxSetF \times supIdxSetC to create supIdxVect
        supIdxVect = reshape(repmat(supIdxVectF,mDimC,1),nDim,mDim);
        supIdxVectC = reshape(repmat(supIdxVectC,1,mDimF),nDim,mDim); 
        supIdxVect = supIdxVect + supIdxVectC;
        supIdxVect(binaryIdx,:) = spones(supIdxVect(binaryIdx,:)); 
        supIdxVect = unique(supIdxVect','rows','first','legacy')'; 
    end
    supIdxVect = grevlexSort(supIdxVect')';
end