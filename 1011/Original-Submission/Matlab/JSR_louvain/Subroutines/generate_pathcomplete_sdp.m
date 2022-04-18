function [A,Agamma,b,c,K] = generate_pathcomplete_sdp(M,graph,options)
%
% GENERATE_PATHCOMPLETE_SDP         returns the SDP problem required for
%                                   finding a common picewise quadratic 
%                                   Lyapunov function for a particular  
%                                   (path-complete) graph.
%
%       NOTE 1 :    The graph must be path-complete to guarantee that the
%                   solution is an upper estimation of the JSR. 
%       NOTE 2 :    This is a low-level function and the validity of inputs 
%                   is not checked.
%
%     [A,Agamma,b,c,K] = GENERATE_PATHCOMPLETE_SDP(M,graph)
%       returns matrices and vectors which describes the following SDP :
%
%             max_P 0
%       s.t.    gamma^2( Mk' Pj Mk ) - Pi<=0 forall edges i,j,k (SDP sense)
%               Pi                      >  0 forall nodes i     (SP sense)
%
%       To describe the above problem for a given gamma with the given
%       output, we have to write :
%
%             min_x c'x
%       s.t.    [ A + (gamma^2-1)Agamma ] x = b
%               x in K
%
%       To have more information about the structure K, type HELP SEDUMI.
%       M must be a cell array of matrices and graph must be a structure
%       with the following fields :
%           graph.nNodes    must contain the number of nodes in the graph.
%           graph.edges     is a nedges x 3 matrix, where the edge i,j,k
%                           describes the following inequality :
%                           gamma^2 ( Ak' Pi Ak ) - Pj <= 0 (SPD sense).
%
%
%     [A,Agamma,b,c,K] = GENERATE_PATHCOMPLETE_SDP(M,graph,options)
%       Does the same as above but with specified parameters described in
%       fields of the structure options.pathcomplete. See JSRSETTINGS and 
%       below for available parameters and options.
%
%       No options are available for the moment.
%
% 
%     WARNING : in [1], we need to have a path-complete graph to have
%       the guarantee that we get an upper bound for the JSR. Other types
%       of graph can be used but without any guarantee.
%
% REFERENCES
%   [1] Ahmadi, Jungers, Parrilo and Roozbehani,
%   "Joint spectral radius and path-complete graph Lyapunov functions"
%   Vol. 52, No1, pp. 687-717 in SIAM J. CONTROL OPTIM, 2014.



matrixSize = size(M{1},1);
nEdges = size(graph.edges,1);
nNodes = graph.nNodes;

if iscellcomplex(M)
    error('Works only with real matrix.')
end

% Index map : symmetric matrix to vector. mapIdx(i,j) gives the index in a
% vector of the (i,j)th element.
% Example : for 2x2 matrices;
%
%   map = [ 1 2 ;
%           2 3 ] ;
%
% It means that, is we have a vector v = [a,b,c], then
% A(1,1) = v(map(1,1)) = v(1) = a ;
% A(2,1) = v(map(2,1)) = v(2) = b ;
% A(1,2) = v(map(1,2)) = v(2) = b ;
% A(2,2) = v(map(2,2)) = v(3) = c ;
%
% Note that for SDP solver, the cone S is for symetric semi-definite
% matrix, that's why the matrix map is also symetric and there is only
% n(n+1)/2 elements in A and not n^2.


mapIdx = zeros(matrixSize);
counter = 0;
for i=1:matrixSize;
    counter = counter + 1;
    mapIdx(i,i) = counter;
    for j=(i+1):matrixSize; % With SDP solver, we assume that P(i,j) = P(j,i)
        counter = counter + 1;
        mapIdx(i,j) = counter;
        mapIdx(j,i) = counter;
    end
end

% Description of the dual problem -> max b'y : AY-C in K* 
% <=> max b'y : AY-C = S, S in K*.
nFreeVarPrimal =  0;
nPositiveVarPrimal = 0;
nSdpMatricesPrimal = nNodes;
nConstraintsPrimal = (nEdges + nNodes) * matrixSize^2;
nVarPrimal = nFreeVarPrimal + nPositiveVarPrimal + nSdpMatricesPrimal*matrixSize*(matrixSize+1)/2;% Symetric!


% Description of the primal problem ( min c'x : Ax = b, X in K )
nFreeVarDual =  0;
nPositiveVarDual = 0;
nSdpMatricesDual = nEdges + nNodes;

% Description of the cone K (primal), see SEDUMI
K.f = nFreeVarDual;
K.l = nPositiveVarDual;
K.q = 0;
K.r = 0;
K.s = ones(nSdpMatricesDual,1)*matrixSize;

% non-zero element in the matrix.
nnz_max = nEdges*(matrixSize^4+matrixSize^2) + nNodes*(matrixSize^2);
nnz_max_gamma = matrixSize^4*nEdges;  % Idem with Agamma

% Check memory used
memUsed = nConstraintsPrimal + nVarPrimal + (nnz_max + nnz_max_gamma)*3; % b + c + A + Agamma
memUsed = memUsed*8; % suppose double.
memNeeded = 0.85*available_memory; % security factor.
if(memUsed > memNeeded)
    error(['Too much memory needed. Memory needed : [',num2str(round(memUsed/1e6)),'] mB. Memory available : ',num2str(round(memNeeded/1e6)) , ' mB.'] )
end
if(memUsed > 0.5*memNeeded)
    warning('The memory used is more than (1/2)* available memory. The program continues but the solver can potentially use more memory and cause a crash.')
end



% Memory allocation
c = zeros(nConstraintsPrimal,1);
b = zeros(nVarPrimal,1);

 % We store A in a sparse structure.
sparseLine = zeros(nnz_max,1); 
sparseCol = zeros(nnz_max,1);
sparseEntry = zeros(nnz_max,1);

sparseLineGamma = zeros(nnz_max_gamma,1);
sparseColGamma = zeros(nnz_max_gamma,1);
sparseEntryGamma = zeros(nnz_max_gamma,1);


% Initialization of some indices
offsetSdpIdx = nFreeVarPrimal + nPositiveVarPrimal ;
constraintIdx = 0;
sparseIdx = 0;
sparseGammaIdx = 0;

% First loop : We describe, for a given label M, M' Pj M - Pi <= 0.
% It can be shown that it is equivalent to
% [ sum[m,n] ( M*_{l,m} M_{n,k} Pj_{m,n} ) - Pi{k,l} ] <= 0 (sdp sense)
%
% Do not forget that P is symetric.

for edgeIdx = 1:nEdges
    i = graph.edges(edgeIdx,1);
    j = graph.edges(edgeIdx,2);
    label = graph.edges(edgeIdx,3);
    Mstar = M{label}';
    
    % for each elements in Pi
    for k=1:matrixSize
        for l=1:matrixSize % We do not need to care about symetry here,
                           % everything is done by the mapIdx matrix.
                                     
            constraintIdx = constraintIdx + 1;
            
            % for each elements in Pj
            for m = 1:matrixSize
                for n = 1:matrixSize % We don't need to care about symetry here,
                                     % everything is done by the mapIdx matrix.
                    
                    Acol = getACol (j,m,n,matrixSize,offsetSdpIdx); % index of variable A(i,j)
                    
                    sparseIdx = sparseIdx+1;
                    sparseLine(sparseIdx) = constraintIdx;
                    sparseCol(sparseIdx) = Acol;
                    sparseEntry(sparseIdx) = Mstar(l,m) * ...% conjugate
                        M{label}(n,k); 
                    
                    % for gamma block
                    sparseGammaIdx = sparseGammaIdx+1;
                    sparseLineGamma(sparseGammaIdx) = constraintIdx;
                    sparseColGamma(sparseGammaIdx) = Acol;
                    sparseEntryGamma(sparseGammaIdx) = sparseEntry(sparseIdx);
                end
            end
            Acol = getACol (i,k,l,matrixSize,offsetSdpIdx);
            
            sparseIdx = sparseIdx+1;
            sparseLine(sparseIdx) = constraintIdx;
            sparseCol(sparseIdx) = Acol;
            sparseEntry(sparseIdx) = -1;
        end
    end
end

% Second loop : we want to describe P > 0. However, we have only a SDP solver.
% So, we use the homogeneity property of the problem (if P is a
% solution, then alpha*P is also a solution, for any positive alpha). It 
% implies that writing P >= I is equivalent to write P > 0.
% NOTE : we have to write the constraint in the dual form - so we  
% write -P - (- I) <= 0, where C = -I.

for i = 1:nNodes
    for k = 1:n
        for l = 1:n
            constraintIdx = constraintIdx+1;
            Acol = getACol (i,k,l,matrixSize,offsetSdpIdx);
            
            sparseIdx = sparseIdx+1;
            sparseLine(sparseIdx) = constraintIdx;
            sparseCol(sparseIdx) = Acol;
            sparseEntry(sparseIdx) = -1;
            if k==l
                c(constraintIdx) = -1;
            end
        end
    end
end

A = sparse(sparseLine,sparseCol,sparseEntry,nConstraintsPrimal,nVarPrimal,nnz_max);
Agamma = sparse(sparseLineGamma,sparseColGamma,sparseEntryGamma,nConstraintsPrimal,nVarPrimal,nnz_max_gamma);

% Col and line are reversed because it is A in the primal form, not A'.
A = A';
Agamma = Agamma';

    function out = getACol (numA,i,j,n,offset)
        nElemInSymMatrix = n*(n+1)/2;
        out = offset+(numA-1)*nElemInSymMatrix + mapIdx(i,j);
    end
end
