function [Q0vec, H0vec, Qpvecs, PSDcone, polyCone, UbdIX] ...
    = BBPOPtoCOP(objPoly, ineqPolySys, I01, Icomp, params)
%BBPOPTOCOP Generate a DNN relaxation problem of the following problem
%
%      zeta* = min. evalPoly(objPoly, x)
%      s.t. evalPoly(ineqPolySys{ii}, x) = 0  (for all ii)
%           x_i \in [0,1]   ( i \in Ibox )
%           x_j \in {0,1}   ( j \in I01 )
%           x^\alpha = 0    ( \alpha = Icomp(ii,:) )
% Usage:
%   [Q0vec, H0vec, Qpvecs, PSDcone, polyCone, UbdIX] = ...
%       BBPOPTOCOP(objPoly, ineqPolySys, I01, Icomp, params);
%
% Note:
%   The input ineqPolySys (which must be a SOS function)
%   and the corresponding output Qpvecs are used 
%   for the Lagrangian-DNN relaxation.


%% Check Inputs
% arrange
numVar = size(objPoly.supports, 2);
I01 = logical(I01);
Icomp = logical(Icomp);
ineqPolySys = ineqPolySys(:);
if isempty(I01); I01 = logical(sparse(0, numVar)); end
if isempty(Icomp); Icomp = logical(sparse(0, numVar)); end
if isempty(ineqPolySys);  ineqPolySys = cell(0);  end
if nargin==4 %|| isempty(params)
    params = defaultParamsDNN();
else
    params = fillUnspecifiedParams(params, defaultParamsDNN());
end

if params.printyesDNN >= 1
    fprintf('\n----------------------------------------------'); 
    fprintf('\n BBPOPtoCOP');
    fprintf('\n----------------------------------------------\n');
end

% format
if ~checkPoly(objPoly); error('objPoly does not follow the sparsePOP format'); end
if ~iscell(ineqPolySys)
    error('ineqPolySys must be a type of cell.');
else
    for ii = 1:length(ineqPolySys)
        if ~checkPoly(ineqPolySys{ii})
            error(['ineqPolySys{' num2str(ii) '} is not correct']);
        end
    end
end

% dim
if length(I01) ~= numVar;  error('I01: Dimensions must agree.'); end
if size(Icomp, 2) ~= numVar;  error('Icomp: Dimensions must agree.'); end
for ii = 1:length(ineqPolySys)
    if size(ineqPolySys{ii}.supports, 2) ~= numVar
        error(['ineqPolySys{' num2str(ii) '}: Dimensions must agree.']);
    end
end

% simplify (reduce degree of binary variables, gather same monomials)
objPoly = simplifyPoly(objPoly, I01, false);
for ii = 1:length(ineqPolySys)
    ineqPolySys{ii} = simplifyPoly(ineqPolySys{ii}, I01, false);
end

% order
if params.relaxOrder < (0.5 * max(sum(objPoly.supports, 2)))
    error('params.relaxOrder is too small to represent objective function.');
end
if params.relaxOrder < (0.5 * max(sum(Icomp, 2)))
    error('params.relaxOrder is too small to represent complementarity.')
end
for ii = 1:length(ineqPolySys)
    if params.relaxOrder < (0.5 * max(sum(ineqPolySys{ii}.supports, 2)))
        error('params.relaxOrder is too small to represent complementarity.')
    end
end

%% Generate Outputs

%ToDo?: [ineqPolySys, I01, Icomp] = reduceConstr(ineqPolySys, I01, Icomp);

if params.printyesDNN >= 2
    fprintf('Generating polyCone...');
    stime = tic;
end
[polyCone, ineqPolySys] = genPolyCone(objPoly, ineqPolySys, I01, Icomp, params); % \mathbb{L}
if params.printyesDNN >= 2
    etime = toc(stime);
    fprintf('Done. (%1.2e sec.)\n', etime);
end

if params.printyesDNN >= 2
    fprintf('Generating PSDcone...');
    stime = tic;
end
PSDcone = genPSDcone(polyCone.varStructure);
if params.printyesDNN >= 2
    etime = toc(stime);
    fprintf('Done. (%1.2e sec.)\n', etime);
end

if params.printyesDNN >= 2
    fprintf('Generating coefficient matrices H0, Q0, Qp...');
    stime = tic;
end
H0vec = genH0vec(polyCone.varStructure);%sparse(1,1,1,length(Q0vec),1);
Q0vec = genQ0vec(polyCone.varStructure, objPoly);
[~, Qpvecs] = genH1vec(polyCone.varStructure, ineqPolySys, params);
if params.printyesDNN >= 2
    etime = toc(stime);
    fprintf('Done. (%1.2e sec.)\n', etime);
end

if params.printyesDNN >= 2
    fprintf('Computing UbdIX...\n');
    stime = tic;
end

if isfield(params,'UbdIX')
    UbdIX = params.UbdIX;
else
    UbdIX = min(estimateUbdIX(polyCone, Icomp, params),...
        estimateUbdIXpd(polyCone, Icomp, params));
end
if params.printyesDNN >= 2
    etime = toc(stime);
    fprintf('Done. (%1.2e sec.)\n', etime);
end

end

%%
function [H1vec, Qpvecs] = genH1vec(varStructure, ineqPolySys, params)
    %idxSupports = varStructure.idxSupports;
    %sizeblk = cellfun(@(x) size(x, 1), idxSupports);
    BMat = varStructure.BMat;
    sumBMat = sum(BMat, 1);
    mm = length(ineqPolySys);
    Qpvecs = zeros(size(BMat, 1), mm);
    %Qpvecs = cell(1, mm);
    for ii = 1:mm
        avecoef = ineqPolySys{ii}.coef ./ sumBMat(varStructure.sysic{ii})';
        Qpvecs(:, ii) = BMat(:, varStructure.sysic{ii}) * avecoef;
        %Qpvecs{ii} = BMat(:, varStructure.sysic{ii}) * sparse(avecoef);
    end
    if mm == 0
        H1vec = sparse(size(BMat, 1), 1);
    else
        H1vec = sum(Qpvecs, 2);
    end
end

%%
function UbdIX = estimateUbdIX(polyCone, Icomp, params)
    sizeblk = polyCone.varStructure.sizeblk;
    if params.LagDNN
        UbdIX = sum((0.5 * (sizeblk - 1)) + 1);
        return
    end
    
    eta = max(sum(Icomp, 1));

    vS = polyCone.varStructure;
    isDiag = vS.rowIdx == vS.colIdx;
    numAppearInDiag = accumarray(vS.momentic(isDiag), 1, [length(vS.ia), 1]);
    cardiB = sum(numAppearInDiag .* ~polyCone.eq0Var);
    if eta == 0
        UbdIX = cardiB; %sum(numAppearInDiag .* ~eq0Var);
        return
    end
    
    
    isNNZinDiag = logical(numAppearInDiag) & ~polyCone.eq0Var;
    numAppearNNZDiag = numAppearInDiag(isNNZinDiag);
    diagSupports = vS.supports(isNNZinDiag, :);
    c = @(S) sum(any(diagSupports(:, S), 2) .* numAppearNNZDiag);
    f = @(S) sum(any(Icomp(:, S), 2));
    
    [~, val] = GSC(c, f, size(Icomp, 2));
    val = full(val);
    
    Heta = sum(1 ./ (1:eta));
    dd = full(max(sum(spones(diagSupports), 2)));
    LB = ceil(val / (dd * Heta));
    UbdIX = floor(cardiB - LB);
    if params.printyesDNN >= 2
        fprintf('cardiB: %d, UbdIX:%d, LB:%d, val:%d, dd:%d, Heta:%g\n', cardiB, UbdIX, LB, val, dd, Heta);
    end
end

%%
function [newObjPoly, newIneqPolySys, newIdxSupports, newI01, newIcomp, xIdx] ...
    = addSlack(objPoly, ineqPolySys, idxSupports, I01, Icomp)
    %ADDSLACK   Add slack variable 
    %

    numXVar = length(I01);

    [uniqueIdx, ~, icIdx] = fastUnique(cell2mat(idxSupports));
    originalBlkSize = cellfun(@(x) size(x, 1), idxSupports);
    icIdx = mat2cell(icIdx, originalBlkSize, 1);
    numSlk = size(uniqueIdx, 1) - 1;
    numSlkInBlk = cellfun(@(x) length(x), icIdx) - 1;

    %% poly
    newObjPoly.coef = objPoly.coef;
    newObjPoly.supports = [objPoly.supports, sparse(size(objPoly.supports, 1), numSlk)];
    
    mm = length(ineqPolySys);
    newIneqPolySys = cell(mm + numSlk, 1);
    for ii = 1:mm
        newIneqPolySys{ii}.coef = ineqPolySys{ii}.coef;
        newIneqPolySys{ii}.supports = ...
            [ineqPolySys{ii}.supports, sparse(size(ineqPolySys{ii}.supports, 1), numSlk)];
    end
    %% addeq
    % is it ok?
    for ii = 1:numSlk
        newIneqPolySys{ii + mm}.coef = [-1; 1; 1];
        newIneqPolySys{ii + mm}.supports = [sparse(1, numXVar + numSlk); %const
            uniqueIdx(ii + 1, :), sparse(1, numSlk); % x^alpha
            sparse(1, numXVar + ii, 1, 1, numXVar + numSlk)]; % slk_alpha
    end

    %% idxSupports
    newIdxSupports = cell(length(idxSupports), 1);
    for ii = 1:length(icIdx)
        newIdxSupports{ii} = ...
            [idxSupports{ii}, sparse(originalBlkSize(ii), numSlk); ...
            sparse(numSlkInBlk(ii), numXVar), ...
            sparse(1:numSlkInBlk(ii), icIdx{ii}(2:end) - 1, 1, numSlkInBlk(ii), numSlk)];
    end
    %newBlkSize = cellfun(@(x) size(x, 1), newIdxSupports);

    %% I01
    slkI01 = ~any(uniqueIdx(2:end, ~I01), 2)'; % binary = there is no non-binary variable
    newI01 = [I01, slkI01];

    %% Icomp
    forcomp = find(slkI01);
    newIcomp = logical([Icomp, sparse(size(Icomp, 1), numSlk);...
        uniqueIdx(forcomp + 1, :),...
        sparse(1:length(forcomp), forcomp, 1, length(forcomp), numSlk)]);

    %% xIdx
    xIdx = [true(1, numXVar), false(1, numSlk)];
end

%%
function PSDcone = genPSDcone(varStructure)
    PSDcone = varStructure.sizeblk;
end

%%
function [polyCone, ineqPolySys] = genPolyCone(objPoly, ineqPolySys, I01, Icomp, params)
    numVarX = size(objPoly.supports, 2);
    eqSupports = cellfun(@(poly) poly.supports, ineqPolySys, 'UniformOutput', false);
    allSupports = [objPoly.supports; Icomp; cell2mat(eqSupports)];
    cliques = genCliques(allSupports, params); % mathcal{C}_k
    
%     if params.sparseSW && params.cliqueMIP
%         cliques = cliqueMIP(cliques, 100);
%     end
    
    sparseIdxSupports = genIdxSupports(cliques, params.relaxOrder, I01); % mathcal{A}_k
    sparseEigTime = sum(cellfun(@(x) size(x, 1), sparseIdxSupports).^3);

    denseIdxSupports = genIdxSupports(true(1, numVarX), params.relaxOrder, I01);
    denseEigTime = size(denseIdxSupports{1}, 1)^3;
    if denseEigTime < sparseEigTime
        idxSupports = denseIdxSupports;
    else
        idxSupports = sparseIdxSupports;
    end

    if params.LagDNN 
        [objPoly, ineqPolySys, idxSupports, I01, Icomp, xIdx] ...
            = addSlack(objPoly, ineqPolySys, idxSupports, I01, Icomp);
    else
        xIdx = true(1, numVarX);
    end
    
    lowTriMomentMat = genLowTriMomentMat(idxSupports, I01); % blk,row,col,support
    varStructure = findCommonMonos(lowTriMomentMat, objPoly, ineqPolySys, numVarX);
    varStructure.sizeblk = countSizeblk(idxSupports);
    varStructure.BMat = genBasisMat(varStructure);
    varStructure.idxSupports = idxSupports;
    varStructure.xIdx = xIdx;
    
    [chain,eq0Var,lbd,ubd] ...
        = genEqIneq(varStructure.supports, I01, Icomp, params.DNNSW);
    if params.findArborescence
        chain = findLongChain(varStructure.supports, eq0Var, chain);
        chain = findArborescence(varStructure.supports, eq0Var, chain);
    elseif params.findLongChain
        chain = findLongChain(varStructure.supports, eq0Var, chain);
    end
    polyCone.chain = chain;
    polyCone.singleton = ~any(chain, 2);
    polyCone.eq0Var = eq0Var;
    polyCone.lbd = lbd; 
    polyCone.ubd = ubd; %???
    polyCone.nonneg = params.DNNSW; %for sdp solver
    polyCone.varStructure = varStructure;
    polyCone.LagDNN = params.LagDNN;
    
    % for proj
    [permIdx, ~] = find(chain);
    cardinality = full(sum(varStructure.BMat, 1))';
    polyCone.varStructure.cardinality = cardinality;
    polyCone.cardinalityPerm = cardinality(permIdx)';

    polyCone.BMatPerm = varStructure.BMat(:, permIdx);
    
    lengthEachChain = full(sum(chain, 1));
    polyCone.lengthEachChain = lengthEachChain;
    polyCone.I = cumsum([0, lengthEachChain]);
    polyCone.topIdx = polyCone.I + 1;
    polyCone.tailIdx = lengthEachChain;
    polyCone.permIdx = permIdx;
end

%%
function sizeblk = countSizeblk(idxSupports)
    numblk = length(idxSupports);
    sizeblk = zeros(1, numblk);
    for ii = 1:numblk
        sizeblk(ii) = size(idxSupports{ii}, 1);
    end
end

%%
function varStructure = findCommonMonos(momentMat, objPoly, ineqPolySys, numVarX)

    numVar = size(objPoly.supports, 2); % with slack
    eqSupports = cellfun(@(poly) poly.supports, ineqPolySys, 'UniformOutput', false);
    [varStructure.supports, varStructure.ia, ic] =...
        fastUnique([momentMat.supports; objPoly.supports; cell2mat(eqSupports); ...
        speye(numVarX), sparse(numVarX, numVar - numVarX)]);

    numMoment = size(momentMat.supports, 1);
    numObj = size(objPoly.supports, 1);
    numSys = cellfun(@(supp) size(supp, 1), eqSupports, 'UniformOutput', false);
    
    numBlk = [numMoment; numObj; cell2mat(numSys); numVarX];
    icBlk = mat2cell(ic, numBlk, 1);

    %varStructure.momentic = ic(1:numMoment);
    %varStructure.objic = ic(numMoment + (1:numObj));
    %varStructure.xic = ic(numMoment + numObj + (1:numVarX));
    varStructure.momentic = icBlk{1};
    varStructure.objic = icBlk{2};
    varStructure.sysic = cell(length(numSys), 1);
    for ii = 1:length(numSys)
        varStructure.sysic{ii} = icBlk{ii+2};
    end
    varStructure.xic = icBlk{end};
    varStructure.rowIdx = momentMat.rowIdx;
    varStructure.colIdx = momentMat.colIdx;
    varStructure.blkIdx = momentMat.blkIdx;
    varStructure.computeMomentIdx = momentMat.computeMomentIdx;

end

%%
function Q0vec = genQ0vec(varStructure, objPoly)
    sumBMat = sum(varStructure.BMat, 1);
    avecoef = objPoly.coef ./ sumBMat(varStructure.objic)';
    Q0vec = varStructure.BMat(:, varStructure.objic) * sparse(avecoef);
end

%%
%
function H0vec = genH0vec(varStructure)
    bb = varStructure.BMat(:, 1);
    H0vec = (1 / sum(bb)) * bb;
end