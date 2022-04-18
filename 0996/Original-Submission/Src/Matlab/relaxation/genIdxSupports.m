function idxSupports = genIdxSupports( cliques, relaxOrder, Ibin )
%GENIDXSUPPORTS   Generate indices of moment matrix.
%
%   This function generates the indices of moment matrix, 
%   (i.e., \mathcal{A}^k_omega)
%
%   INPUT:
%      cliques    : the k-th row is a logical vector that specifies the
%                   variables in the k-th clique (i.e., V^k)
%      relaxOrder : relaxation order of DNN hierarchy.
%      Ibin       : logical vector which specifies binary variables.
%
%   OUTPUT:
%      idxSupports: type of cell. k-th cell contains two dimensional array
%                   whose row represents the elements of \mathcal{A}_k.
%
% Reference:
% N. Ito, S. Kim, M. Kojima, A. Takeda, and K.-C. Toh.
% A Sparse Doubly Nonnegative Relaxation of Polynomial Optimization Problems
% with Binary, Box and Complementarity Constraints,
% Research Rport B-48?, Department of Mathematical and Computing Sciences, 
% Oh-Okayama, Meguro-ku, Tokyo 152-8552, March 2018. 


numVar = size(cliques,2);

idxSupports = cell(size(cliques,1),1);
for ii = 1:size(cliques,1)
    iiClique = logical(cliques(ii,:));

    iiBoxVar = iiClique & ~Ibin;
    iiBinVar = iiClique & Ibin;
    numiiBoxVar = sum(iiBoxVar);
    numiiBinVar = sum(iiBinVar);
    
    tmp = enumMonosUptoDeg(numiiBoxVar+1, relaxOrder);
    trashIdx = tmp(:,1) > numiiBinVar;
    degBin = tmp(~trashIdx, 1);
    boxMonos = tmp(~trashIdx, 2:end);
    
    binMonos = enum01MonosUptoDeg(numiiBinVar,max(degBin));
    sumBin = sum(binMonos,2);
    
    ddmax = max(degBin);
    boxIdx = cell(ddmax+1,1);
    binIdx = cell(ddmax+1,1);
    for dd = 0:ddmax
        vec1 = find(degBin == dd);
        vec2 = find(sumBin == dd);
        [ddBoxIdx,ddBinIdx] = myMeshGrid(vec1,vec2);
        boxIdx{dd+1} = ddBoxIdx;
        binIdx{dd+1} = ddBinIdx;
    end
    boxIdx = cell2mat(boxIdx);
    binIdx = cell2mat(binIdx);
    
    iiSupports = sparse(length(boxIdx),numVar);
    iiSupports(:,iiBoxVar) = boxMonos(boxIdx,:);
    iiSupports(:,iiBinVar) = binMonos(binIdx,:);
    idxSupports{ii} = grevlexSort(iiSupports);
end
end

% faster than meshgrid
function [Row,Col] = myMeshGrid(vec1,vec2)
    if ~isrow(vec1); vec1 = vec1'; end
    if ~isrow(vec2); vec2 = vec2'; end
    
    Row = ones(length(vec2),1) * vec1;  % faster than repmat
    Col = vec2' * ones(1,length(vec1)); % faster than repmat
    Row = Row(:);
    Col = Col(:);
end