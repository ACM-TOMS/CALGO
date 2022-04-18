function cliques = genCliques(supports, params)
%GENCLIQUES   Generate cliques in correlation sparsity matrix
%   Usage:
%      cliques = genCliques(supports, params)
%   Input:
%      supports : supports of polynomials in sparsePOP format
%      params.sparseSW : logical value to specify if it exploits sparsity
%      params.cliques  : you can manually specify cliques
%   Output:
%      cliques: boolean matrix which represents maximal cliques.
%         Each row corresponds to each maximal clique.
%         Each column corresponds to variable (vertex).
%         For example, if there are vertices {1,2,3} and
%         C1={1,2} and C2={2,3} are maximal cliques, then
%         cliques = logical([1 1 0; 0 1 1]);

    if ~params.sparseSW
        cliques = 1:size(supports, 2);
        return;
    end
    
    if isfield(params,'cliques')
        %checkCliques(params.cliques)
        cliques = params.cliques;
        return;
    end
    
    numVar = size(supports, 2);
    
    %% Step 1:
    % Generate csp matrix by gathering all the supports.
    uniqueSupports = fastUnique(spones(supports));
    numTerms = size(uniqueSupports, 1);
  
    tuniqueSupports = logical(uniqueSupports)';
    
    rowIdx = cell(numTerms, 1);
    colIdx = cell(numTerms, 1);
    for ii = 1:numTerms
      nnzIdx = find(tuniqueSupports(:, ii));
      [rowIdx{ii}, colIdx{ii}] = myMeshGrid(nnzIdx, nnzIdx);
    end
    rowIdx = cell2mat(rowIdx);
    colIdx = cell2mat(colIdx);
    cspmat = spones(sparse(rowIdx, colIdx, 1, numVar, numVar));
    
    % to positive definite (diagonally dominant)
    cspmat = cspmat + sparse(1:numVar, 1:numVar, sum(cspmat, 2) + 0.1); 
  
  
    %% Step 2
    % Chordal extension by Cholesky decomposition
  
    I = symamd(cspmat); % Pre-process: Minimum degree ordering
    [R, p] = chol(cspmat(I, I));
    if (p > 0) 
      error('Correlative sparsity matrix is not positive definite.');
    end
    
    %% Step3
    % Find the maxmal clieques in the chodal graph
    R = logical(R);
    [~, revI] = sort(I);
    isNewMaxClique = false(numVar, 1);
    isNewMaxClique(1) = true;
    for ii = 2:numVar
      isNotMaximal = any(all(R(1:(ii - 1), R(ii, :)), 2));
      isNewMaxClique(ii) = ~isNotMaximal;
    end
    cliques = R(isNewMaxClique, revI);
end