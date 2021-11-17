function [objPoly,I01,Icomp,relaxOrder,params] = ...
    genPOPdense(degree,nDim,isBin,addComplement)
%% [objPoly,I01,Icomp,relaxOrder,params] = genPOPdense(degree,nDim,isBin,addComplement)
%% Input
%%   degree  : the degree of the polynomial objective function
%%   nDim    : the dimension of the polynomial objective function
%%   isBin   : 1 for the binary constraint, 0 for the box constraint
%%   addComplement   : 1 for adding random complementarity constraints 

    relaxOrder = ceil(degree/2); 
    params.delta1 = 1.0e-4;
    params.delta = 0;
    
    randSeed = 2017;
    I01 = isBin .* true(1, nDim);
    probData.nDim = nDim;
    probData.degree = degree;
    probData.randSeed = randSeed;
    probData.density = 1;
    %%%%%%%%%%
    rng('default');
    rng(probData.randSeed);
    objPoly.typeCone = 1;
    objPoly.sizeCone = 1;
    objPoly.dimVar = probData.nDim;
    objPoly.degree = probData.degree;
    % enumerate all supports
    objPoly.supports = ...
        fullSup2(probData.nDim, find(I01), find(~I01), objPoly.degree)';
    if isfield(probData,'density') && probData.density < 1
        % remove supports
        dVect = sum(objPoly.supports,2); 
        idxAll = [];
        for degree=2:objPoly.degree 
            idxEach = find(dVect == degree);
            rVect = spones(sprand(1,length(idxEach),probData.density)); 
            I = find(rVect); 
            if ~isempty(I)
                idxEach = idxEach(I);
            end
            idxAll = [idxAll;idxEach(:)]; 
        end
        objPoly.supports = objPoly.supports(idxAll,:); 
        % objPoly.degree = max(dVect); 
    end
    % generate coef
    objPoly.noTerms = size(objPoly.supports,1);
    objPoly.coef = 2*rand(objPoly.noTerms,1) - ones(objPoly.noTerms,1);
    % compute UbdObjVal
    x0 = zeros(objPoly.dimVar,1);
    objPoly.UbdObjVal = evalPoly(objPoly,x0);
    
    clique.Set = true(1, probData.nDim);
    Icomp = [];
    if addComplement
        Icomp = complementarityConstraint(2, 2, clique, randSeed);
    end
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