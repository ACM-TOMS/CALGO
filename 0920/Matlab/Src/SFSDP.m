function [xMatrix,info,distanceMatrix] = SFSDP(sDim,noOfSensors,noOfAnchors,xMatrix0,distanceMatrix,pars)
% 
% All the functions, except sedumi and sedumiwrap (a MATLAB version of SDPA), 
% used in this program are included in the file SFSDP.m
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SFSDP --- A sparse version of FSDP
% Sunyoung Kim, Masakazu Kojima^* and Hayato Waki
% July , 2008
%
% Revised July 2009
% Sunyoung Kim, Masakazu Kojima^*, Hayato Waki and Makoto Yamashita
%
% * Department of Computing and Mathematical Sciences
%   Tokyo Institute of Technology
%   Oh-Okayama, Meguro, Tokyo 152-8552
%   e-mail: kojima@is.titech.ac.jp
%
% SFSDP is proposed in the paper
% S. Kim, M. Kojima and H. Waki, "Exploiting Sparsity in SDP Relaxation 
% for Sensor Network Localization", SIAM Journal of Optimization
% Vol.20 (1) 192-215 (2209). 
%
% SFSDP is a sparse version of FSDP which was proposed in 
% P. Biswas and Y. Ye (2004) ÅgSemidefinite programming for ad hoc wireless 
% sensor network localization,Åh in Proceedings of the third international 
% symposium on information processing in sensor networks, ACM press, 46-54. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Acknowledgment
% 
% The authors of this software package are grateful to Professor Yinyu Ye
% who provided us with the original version of FSDP, and Professor Kim Chuan Toh 
% for lots of helpful comments and providing MATLAB programs refineposition.m 
% and procrustes.m. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is a component of the software package SFSDP 
% Copyright (C) 2008 Sunyoung Kim, Masakazu Kojima and Hayato Waki
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sDim :    the dimension of the space in which sensors and anchors are
%           located; sDim is either 2 or 3. 
% noOfSensors   : the number of sensors.
% noOfAnchors   :   the number m_a of anchors; the last m_a columns of xMatrix0, 
% xMatrix0  :   sDim x n matrix of sensors' and anchors' locations in the sDim dimensional space, 
%               where n is the total number of sensors and anchors, and 
%               anchors are placed in the last m_a columns,                
%           or
%               sDim x m_a matrix of anchors in the sDim dimensional space,
%               where m_a denotes the number of anchors. 
%               If noOfAnchors == 0 then xMatrix0 can be [];
% distanceMatrix : the sparse (and noisy) distance matrix between xMatrix0(:,i) and xMatrix0(:,j);
%               distanceMatrix(i,j) = the (noisy) distance between xMatrix0(:,i) and xMatrix0(:,j) if i < j; 
%               distanceMatrix(i,j) = 0 if i >= j.
% pars --- parameters
%   parameters used in SeDuMi and SDPA, default values: 
%       pars.eps = 1.0e-5 for sedumi and sdpa;
%       pars.free = 0 for sedumi;
%       pars.fid = 0 for sedumi and sdpa;
%       pars.alg, pars.theta, pars.beta, pars.vplot, pars.maxiter for sedumi; 
%       See the manual of SeDuMi. 
%   new parameters used in SFSDP
%       pars.SDPsolver  = 'sedumi' to solver the SDP problem by SeDuMi; default
%                       = 'sdpa' to solver the SDP problem by SDPA
%       pars.minDegree = sDim + 2;
%           This parameter control the minimum degree of each sensors. If
%           we increase pars.minDegree, we expect a stronger relaxation but
%           more larger cputime. 
%           Set pars.minDegree >= 100 if no reduction of the edges is
%           necessary. 
%       pars.objSW 
%           = 0 if #anchors >= sDim+1 and no noise 
%                   ---> solving quadratic equalities with objective function = 0
%           = 1 if #anchors >= sDim+1 and noise 
%                   ---> minimizing the 1-norm error
%           = 2 if #anchors = 0  and no noise
%                   ---> minimizing a regularization term 
%           = 3 if #anchors = 0 and noise
%                   ---> minimizing the 1-norm error + a regularization term 
%       pars.noisyFac 
%           = [] if noisyFac is not specified or unknown 
%           = noisyFac if noisFac is known. This parameter will be use to
%               bound the error epsilon_{ij}^+ and epsilon_{ij}^-. 
%       pars.analyzeData 
%           = 1 then the input data is analyzed ---> default
%           = 0 then no information on the input data
%       pars.moreDistanceSW 
%           = 1 then all edges within the radio range are added before
%               constructing a FSDP relaxation
%           = 0 then the original sensor network problem will be solved --- default
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Additional parameters ---> 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       pars,'edgeRemoveSW
%           some different methods for eliminating edges before constructing an SDP
%           relaxation problems
%           1 --- the orignal method
%           2 --- sorting nodes according to their degrees, smaller degree nodes ---> larger degree nodes
%           3 --- sorting nodes according to their degrees, larger degree nodes ---> smaller degree nodes
%           4 --- using the MATLAB function symrcm
%           5 --- rearranging the nodes randomly
%           10 --- choosing the best one among the above methods --- default
%       pars.sparseSW
%           pars.sparseSW = 1 ---> SFSDP
%           pars.sparseSW = 0 ---> FSDP by Biswas and Ye 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% <--- Additional parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% xMatrix   :   sDim x n matrix of sensors' and anchors' locations computed 
%               in the sDim dimensional space, 
%               where n is the total number of sensors and anchors, and 
%               anchors are placed in the last m_a columns 
% info      :   info from SeDuMi output
%   info.cpusec     : the cpu time in seconds for the solution time.
%   info.iter       : the number of iterations. 
%   info.numerr, info.pinf, info.dinf   : see the manual of SeDuMi.
% distanceMatrix    :   the distance submatrix used to construct an SDP
%                       relaxation problem in SFSDP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Additional parameters ---> 
SDPsolverDefault = 'sdpa';
% SDPsolverDefault = 'sedumi';
sparseSWdefault = 1;
%sparseSWdefault = 0;
% <--- Additional parameters

if ~isfield(pars,'SDPsolver')
    pars.SDPsolver = SDPsolverDefault; 
end
if ~isfield(pars,'sparseSW')
    pars.sparseSW = sparseSWdefault; % ---> SFSDP
%    pars.sparseSW = 0; % ---> FSDP by Biswas and Ye 
end

%%%%%%%%%% 
edgeRemoveSWdefault = 10; 
if ~isfield(pars,'edgeRemoveSW')
    pars.edgeRemoveSW = edgeRemoveSWdefault;
    % different methods for eliminating edges before constructing an SDP
    % relaxation problems
    % 1 --- the orignal dfault method
    % 2 --- sorting nodes according to their degrees, smaller degree nodes ---> larger degree nodes
    % 3 --- sorting nodes according to their degrees, larger degree nodes ---> smaller degree nodes
    % 4 --- using the MATLAB function symrcm
    % 5 --- rearranging the nodes randomly
    % 10 --- choosing the best one among the above methods
end
%%%%%%%%%

if pars.sparseSW == 1
    fprintf('\nSFSDP --- A Sparse version of FSDP (Biswas and Ye)\n');
    fprintf('Sunyoung Kim, Masakazu Kojima, Hayato Waki and Makoto Yamashita\n');
    fprintf('Version 1.11, July 2009\n\n');
elseif pars.sparseSW == 0 
    fprintf('\nA modified version of FSDP (Biswas and Ye)\n');
    fprintf('Sunyoung Kim, Masakazu Kojima, Hayato Waki and Makoto Yamashita\n');
    fprintf('Version 1.11, July 2009\n\n');
end
% startingTime = cputime; 
startingTime = tic; 

% Checking the sizes of xMatrix0 --->
if noOfAnchors > 0
    if isempty(xMatrix0) 
        error('## xMatrix0 needs to be sDim x (noOfSensors + noAnchors) or sDim x noAnchors.');
    elseif size(xMatrix0,1) ~= sDim
        error('## xMatrix0 needs to be sDim x (noOfSensors + noAnchors) or sDim x noAnchors.');
    elseif (size(xMatrix0,2) ~= (noOfSensors + noOfAnchors)) && (size(xMatrix0,2) ~= noOfAnchors) 
        error('## xMatrix0 needs to be sDim x (noOfSensors + noAnchors) or sDim x noAnchors.');
    end
else % noOfAnchors == 0
    if ~isempty(xMatrix0) 
        if (size(xMatrix0,1) ~= sDim) || (size(xMatrix0,2) ~= noOfSensors) 
            error('## xMatrix0 needs to be sDim x (noOfSensors + noAnchors) or sDim x noAnchors.');
        end
    end
end
rmsdSW = 1; 
if noOfAnchors == 0
    if isempty(xMatrix0) 
        fprintf('## an anchor free problem and sensor locations are not given\n'); 
        xMatrix0 = sparse(sDim,noOfSensors); 
        rmsdSW = 0; 
    end
elseif noOfAnchors < sDim+1
    fprintf('## case 1 <= noOfAnchors < sDim+1 can not be handled effectively\n'); 
    error('   take noOfAnchors = 0 and modify distanceMatrix'); 
elseif (noOfAnchors > 0) && (size(xMatrix0,2) == noOfAnchors)
    % Only anchors' location are given ---> expand the distanceMatrix 
	% to the size sDim x (noOfSensors + noOfAnchors)
    rmsdSW = 0; 
    fprintf('## only anchor locations are given\n'); 
    xMatrix0 = [sparse(sDim,noOfSensors), xMatrix0(:,1:noOfAnchors)]; 
end 
% <--- Checking the sizes of xMatrix0
% Checking the size of distanceMatrix ---> 
if (size(distanceMatrix,2) ~= noOfSensors + noOfAnchors) || (size(distanceMatrix,1) ~= noOfSensors)
	error('## distanceMatrix needs to be noAnchors x (noOfSensors + noAnchors)');
end
if norm(tril(distanceMatrix),inf) > 1.0e-8
    fprintf('## distanceMatrix should be upper triangular\n');
    fprintf('   distanceMatrix = triu(distanceMatrix,1)\n');
    distanceMatrix = triu(distanceMatrix,1);
end
% <--- Checking the size of distanceMatrix

if nargin < 6
    pars.eps = 1.0e-5;
    pars.free = 0;
    pars.minDegree = sDim + 2;
    pars.noisyFac = []; 
    if noOfAnchors >= sDim+1
        pars.objSW = 1;
    else
        pars.objSW = 3;
    end
else
    if ~isfield(pars,'eps')
        pars.eps = 1.0e-5;
    end
    if ~isfield(pars,'free')
        pars.free = 0;
    end
    if ~isfield(pars,'minDegree')
        pars.minDegree = sDim+2;
    end
    if ~isfield(pars,'noisyFac')
        pars.noisyFac = [];
    end
    if ~isfield(pars,'objSW')
        if noOfAnchors >= sDim+1
            pars.objSW = 1;
        else
            pars.objSW = 3;
        end
        % pars.eps = 1.0e-5;
    end 
    if ~isfield(pars,'fid')
        pars.fid = 0;
    end
end

fprintf('## pars: eps = %6.2e, free = %d, minDegree = %d, objSW = %d',...
    pars.eps,pars.free,pars.minDegree,pars.objSW);
if (~isempty(pars.noisyFac)) && (pars.noisyFac > 0)
    fprintf(', noisyFac = %6.2e\n',pars.noisyFac);
else
    fprintf('\n'); 
end

noOfAnchors0 = noOfAnchors; 

if noOfAnchors0 == 0
    % Anchor free case --->
    % Fixing 3 or 4 points both in the 2 and 3 dimensional cases, respectively, and
    % adding regularization terms to the objective function later.
    [clique,anchorLocation] = findClique(distanceMatrix,sDim+1);
    controlSW = 0; 
    if ~isempty(clique)
        if isempty(anchorLocation)
            rand('state',2414);
            i=0;
            controlSW = 0; 
            while controlSW == 0
                i = i+1;
                [temp,p1] = sort(rand(1,noOfSensors));
                distanceMatrix2 = distanceMatrix + distanceMatrix';
                distanceMatrix2 = distanceMatrix2(p1,p1);
                distanceMatrix2 = triu(distanceMatrix2,1); 
                [clique,anchorLocation] = findClique(distanceMatrix2,sDim+1);
                if (~isempty(clique)) && (~isempty(anchorLocation))
                    controlSW = 1; 
                elseif i == 10
                    controlSW = -1;
                end
            end
            clear distanceMatrix2
            clique = p1(clique);  
        else
            controlSW = 1;
        end
        if controlSW == 1
            if sDim == 2
                fprintf('## an anchor free problem. Fix 3 sensors location temporarily ---> 3 anchors \n');
            elseif sDim == 3
                fprintf('## an anchor free problem. Fix 4 sensors location temporarily ---> 4 anchors \n');
            end
            distanceMatrix = distanceMatrix + distanceMatrix';
            sensorsIdx = setdiff(1:noOfSensors,clique);
            perm1 = [sensorsIdx,clique];
            [temp,perm0] = sort(perm1);
            distanceMatrix = distanceMatrix(perm1,perm1);
            debugSW = 0;
            if debugSW == 1
                full(xMatrix0)
                full(distanceMatrix)
                xMatrix2 = xMatrix0(:,perm1);
                dMat = sparse(noOfSensors,noOfSensors);
                for i=1:noOfSensors
                    for j=i+1:noOfSensors
                        if distanceMatrix(i,j) > 0
                            dMat(i,j) = norm(xMatrix2(:,i)' -xMatrix2(:,j)');
                        end
                    end
                end
                norm(full(dMat+dMat'-distanceMatrix))
                full(anchorLocation + 0.5)
             end
            distanceMatrix = triu(distanceMatrix,1);
            xMatrix0(:,1:noOfSensors-sDim-1) = sparse(sDim,noOfSensors-sDim-1);
            xMatrix0(:,noOfSensors-sDim:noOfSensors) = anchorLocation + 0.5;
            noOfSensors = noOfSensors-sDim-1;
            noOfAnchors = sDim+1;
            %%%%
            distanceMatDeleted = distanceMatrix(noOfSensors+1:noOfSensors+noOfAnchors,:); 
            %%%%
            distanceMatrix = distanceMatrix(1:noOfSensors,:);
         end
    else
        controlSW = -1;
    end
    if controlSW == -1
        fprintf('## an anchor free problem. Fix 2 sensors location temporarily. \n');
        distanceMatrix = distanceMatrix + distanceMatrix';
        % Fixing two points as anchors --->
        IMat = spones(distanceMatrix + distanceMatrix') + (noOfSensors+1)*speye(noOfSensors,noOfSensors);
        perm = symrcm(IMat);
        ix = perm(fix(noOfSensors/2));
        [dist,iy] = max(distanceMatrix(ix,:));
        % Moving the two sensors to the (noOfSensors-1)th and (noOfSensors)th ancors, and fix their
        % coodinates.
        % the permutation being used to retrieve the original arrangement --->
        perm0 = (1:noOfSensors);
        perm0(ix) = noOfSensors-1;
        perm0(iy) = noOfSensors;
        perm0(noOfSensors-1) = ix;
        perm0(noOfSensors) = iy;
        % <--- the permutation being used to retrieve the original arrangement
        xMatrix0 = xMatrix0(:,perm0);
        distanceMatrix = distanceMatrix(perm0,perm0);
        distanceMatrix = triu(distanceMatrix);
        xMatrix0(:,noOfSensors-1) = 0.5;
        xMatrix0(:,noOfSensors) = 0.5;
        xMatrix0(1,noOfSensors) = xMatrix0(1,noOfSensors) + dist;
        xMatrix0(:,1:noOfSensors-2) = zeros(sDim,noOfSensors-2);
        %xMatrix0(:,1:noOfSensors-2) = 0;
        noOfSensors = noOfSensors-2;
        noOfAnchors = 2;
        %%%%
        distanceMatDeleted = distanceMatrix(noOfSensors+1:noOfSensors+noOfAnchors,:);
        %%%%
        distanceMatrix = distanceMatrix(1:noOfSensors,:);
        % <--- Fixing two points as anchors
    end
% elseif noOfAnchors0 == 1
%     fprintf('## only one anchor is given. fix one sensor location temporarily.\n');
%     fprintf('   the problem is processed as an anchor free case.\n');
%     ix = noOfSensors + 1;
%     %
%     % 2008-06-08 Waki
%     % MATLAB recommend the exploitation of max(distanceMatrix(:,idx),[],2) instead of the
%     % following command:
%     %
%     [dist,iy] = max(distanceMatrix(:,ix)');
%     %
%     % But I have not checked whether two commands are equivalent or not, yet.
%     %
%     if dist == 0
%         error('## the one anchor given is not linked any sensors!');   
%     end
%     % Moving the sensor iy to (noOfSensors)th position, and fix its
%     % coordinates. 
%     % the permutation being used to retrieve the original arrangement --->
%     perm0 = (1:noOfSensors+1); 
%     perm0(iy) = noOfSensors;
%     perm0(noOfSensors) = iy;
%     % <--- the permutation being used to retrieve the original arrangement
%     xMatrix0 = xMatrix0(:,perm0);
%     distanceMatrix = [distanceMatrix;sparse(1,noOfSensors+1)] + [distanceMatrix;sparse(1,noOfSensors+1)]';     
%     distanceMatrix = distanceMatrix(perm0,perm0);
%     distanceMatrix = triu(distanceMatrix); 
%     %%%%
%     distanceMatDeleted = distanceMatrix(noOfSensors+1:noOfSensors+noOfAnchors,:);
%     %%%%
%     distanceMatrix = distanceMatrix(1:noOfSensors,:);
%     xMatrix0(:,noOfSensors) = xMatrix0(:,noOfSensors+1);   
%     xMatrix0(1,noOfSensors) = xMatrix0(1,noOfSensors) + dist; 
%     noOfSensors = noOfSensors-1; 
%     noOfAnchors = noOfAnchors + 1; % = 2; 
%     distanceMatrix = distanceMatrix(1:noOfSensors,:); 
end

% Reduction of the sensor network graph by eliminating (redundant) edges ---> 
% This is a heuristic method to get computational efficiency, no
% theoretical guarantee. 

if (pars.minDegree <= 99) && (pars.edgeRemoveSW == 1)
    % not sorting 
    [distanceMatrix,noOfNonzeros1d,noOfNonzeros1] = ... 
        removeRedundantEdges1(sDim,noOfSensors,distanceMatrix,pars.minDegree); 
elseif (pars.minDegree <= 99) && (pars.edgeRemoveSW == 2)
    % smaller degree nodes ---> larger degree nodes
    [distanceMatrix,noOfNonzeros2d,noOfNonzeros2] = ... 
        removeRedundantEdges2(sDim,noOfSensors,distanceMatrix,pars.minDegree); 
elseif (pars.minDegree <= 99) && (pars.edgeRemoveSW == 3)
    % larger degree nodes ---> smaller degree nodes 
    [distanceMatrix,noOfNonzeros3d,noOfNonzeros3] = ... 
        removeRedundantEdges3(sDim,noOfSensors,distanceMatrix,pars.minDegree); 
elseif (pars.minDegree <= 99) && (pars.edgeRemoveSW == 4)
    % using symrcm
    [distanceMatrix,noOfNonzeros4d,noOfNonzeros4] = ... 
        removeRedundantEdges4(sDim,noOfSensors,distanceMatrix,pars.minDegree); 
elseif (pars.minDegree <= 99) && (pars.edgeRemoveSW == 5)
    % sorting nodes randomly 
    [distanceMatrix,noOfNonzeros5d,noOfNonzeros5] = ... 
        removeRedundantEdges5(sDim,noOfSensors,distanceMatrix,pars.minDegree); 
elseif (pars.minDegree <= 99) && (pars.edgeRemoveSW == 10)
    noOfNz = [];
    [distanceMatrix1,noOfNonzeros1d,noOfNonzeros1,maxCliqueSize1] = ... 
        removeRedundantEdges1(sDim,noOfSensors,distanceMatrix,pars.minDegree); 
    noOfNz = [noOfNz, noOfNonzeros1];
    [distanceMatrix2,noOfNonzeros2d,noOfNonzeros2,maxCliqueSize2] = ... 
        removeRedundantEdges2(sDim,noOfSensors,distanceMatrix,pars.minDegree);     
    noOfNz = [noOfNz, noOfNonzeros2];
    [distanceMatrix3,noOfNonzeros3d,noOfNonzeros3,maxCliqueSize3] = ... 
        removeRedundantEdges3(sDim,noOfSensors,distanceMatrix,pars.minDegree); 
    noOfNz = [noOfNz, noOfNonzeros3];
    [distanceMatrix4,noOfNonzeros4d,noOfNonzeros4,maxCliqueSize4] = ... 
        removeRedundantEdges4(sDim,noOfSensors,distanceMatrix,pars.minDegree); 
    noOfNz = [noOfNz, noOfNonzeros4];
    [distanceMatrix5,noOfNonzeros5d,noOfNonzeros5,maxCliqueSize5] = ... 
        removeRedundantEdges5(sDim,noOfSensors,distanceMatrix,pars.minDegree); 
    noOfNz = [noOfNz, noOfNonzeros5];
    [temp,sortIdx] = sort(noOfNz);
    if sortIdx(1) == 1
        distanceMatrix = distanceMatrix1; 
    elseif sortIdx(1) == 2
        distanceMatrix = distanceMatrix2; 
    elseif sortIdx(1) == 3
        distanceMatrix = distanceMatrix3; 
    elseif sortIdx(1) == 4
        distanceMatrix = distanceMatrix4; 
    elseif sortIdx(1) == 5
        distanceMatrix = distanceMatrix5; 
    end
    debugSW = 0; 
    if debugSW == 1
        fprintf('## noOfNonzeros1d = %d, noOfNonzeros1 = %d, maxCliqueSize1 =%d\n',noOfNonzeros1d,noOfNonzeros1,maxCliqueSize1);
        fprintf('## noOfNonzeros2d = %d, noOfNonzeros2 = %d, maxCliqueSize2 =%d\n',noOfNonzeros2d,noOfNonzeros2,maxCliqueSize2);
        fprintf('## noOfNonzeros3d = %d, noOfNonzeros3 = %d, maxCliqueSize3 =%d\n',noOfNonzeros3d,noOfNonzeros3,maxCliqueSize3);
        fprintf('## noOfNonzeros4d = %d, noOfNonzeros4 = %d, maxCliqueSize4 =%d\n',noOfNonzeros4d,noOfNonzeros4,maxCliqueSize4);
        fprintf('## noOfNonzeros5d = %d, noOfNonzeros5 = %d, maxCliqueSize5 =%d\n',noOfNonzeros5d,noOfNonzeros5,maxCliqueSize5);
        sortIdx
        XXXXX       
    end
    clear distanceMatrix1 distanceMatrix2 distanceMatrix3 distanceMatrix4 distanceMatrix5
    debugSW = 0; 
end

noOfEdges = nnz(distanceMatrix); 
noOfEdges2 = nnz(distanceMatrix(:,1:noOfSensors)); 
degreeVector = sum(spones([distanceMatrix(:,1:noOfSensors)+distanceMatrix(:,1:noOfSensors)',...
    distanceMatrix(:,noOfSensors+1:noOfSensors+noOfAnchors)]),2); 
minDeg = full(min(degreeVector')); 
maxDeg = full(max(degreeVector')); 
averageDeg = full(sum(degreeVector')/noOfSensors); 

fprintf('## the number of dist. eq. used in SFSDP between two sensors  = %d\n',noOfEdges2);
fprintf('## the number of dist. eq. used in SFSDP between a sensor & an anchor = %d\n',noOfEdges-noOfEdges2);
fprintf('## the min., max. and ave. degrees over sensor nodes = %d, %d, %6.2f\n',minDeg,maxDeg,averageDeg);

% 
% <--- Reduction of the sensor network graph by eliminating (redundant) edges

[A,b,c,K] = generatePrimalSDP(xMatrix0,noOfAnchors,distanceMatrix,pars);

% the dimensions and sizes of vectors and matrices used in the FSDP
% relaxations ---> 
mDim1 = nnz(distanceMatrix); % the number of distant equations
% the total number of equalities =
% the row size of the constraint matrix A = the number of equality constraitns
if (pars.objSW == 0) || (pars.objSW == 2)
    nDimLP = 0;
else
    nDimLP = 2*mDim1;
end
% the number of LP variables;
% each nonzero distance yields 2 LP variables

nDimSDP = (noOfSensors+sDim)*(noOfSensors+sDim);
% the size of the vectorized SDP block in the sedumi format
nDim = nDimLP + nDimSDP;
% <--- the dimensions and sizez of vectors and matrices used in the FSDP
% relaxations 
        
% Only for anchor free and cases but not used because not effective --->
% if (noOfAnchors0 == 0) && (sDim == 3)
% if  ((pars.objSW == 2) || (pars.objSW == 3)) && (sDim == 3)
%     % Fixing the 3rd cordinate of another point to 0.5
%     iz = 1;
%     addAMat = sparse(noOfSensors+sDim,noOfSensors+sDim);
% %     addAMat(3,sDim+iz) = 1;
% %     addAMat(sDim+iz,3) = 1;
%     addAMat(noOfSensors+3,iz) = 1;
%     addAMat(iz,noOfSensors+3) = 1;
%     if (pars.objSW == 0) || (pars.objSW == 2)
%         addA = reshape(addAMat,1,nDim);
%     else
%         addA = [sparse(1,nDimLP),reshape(addAMat,1,nDim-nDimLP)];
%     end
%     A = [A; addA];
%     b = [b; 1];
% end
% <--- Only for anchor free cases but not used because not effective

% Adding regularization terms to the objective function for anchor free
% cases --->
if (pars.objSW == 2) || (pars.objSW == 3)
    if (pars.sparseSW == 1)
        % t0 = cputime; 
        perturbeEpsion1 = 0.01;
        perturbeEpsion2 = 1;
        pertbationMat = sparse(noOfSensors+sDim,noOfSensors+sDim);
        rand('state',3201);
        regularizationSW = 2;
        % not used ---> 
        if regularizationSW == 1
            for i=1:noOfSensors
                nzIndices = find(distanceMatrix(i,1:noOfSensors));
                for j =nzIndices
                    r = rand(1,1);
                    pertbationMat(i,i) = pertbationMat(i,i) -perturbeEpsion2*r;
                    pertbationMat(j,j) = pertbationMat(j,j) -perturbeEpsion2*r;
                    pertbationMat(i,j) = pertbationMat(i,j) +perturbeEpsion2*r;
                    pertbationMat(j,i) = pertbationMat(j,i) +perturbeEpsion2*r;
                end
            end
        % <--- not used 
        elseif regularizationSW == 2
            pertbationMat0 = sparse(noOfSensors+sDim,noOfSensors+sDim);
            [sparsityPatternMat] = genSparsityPatternMat(A,c,K);
            sparsityPatternMat = sparsityPatternMat(1:noOfSensors,1:noOfSensors) + (1+noOfSensors)*10*speye(noOfSensors,noOfSensors);
            permutation = symamd(sparsityPatternMat);
            [RMat,p] = chol(sparsityPatternMat(permutation,permutation));
            for i=1:noOfSensors
                 nzIndices = find(RMat(i,1:noOfSensors));
                for j=nzIndices
                    if j > i
                        pertbationMat0(i,i) = pertbationMat0(i,i) -perturbeEpsion2;
                        pertbationMat0(j,j) = pertbationMat0(j,j) -perturbeEpsion2;
                        pertbationMat0(i,j) = pertbationMat0(i,j) +perturbeEpsion2;
                        pertbationMat0(j,i) = pertbationMat0(j,i) +perturbeEpsion2;
                    end
                end
            end
            pertbationMat(permutation,permutation) = pertbationMat0(1:noOfSensors,1:noOfSensors);
            clear pertbationMat0 permutation
        end
        % t1 = cputime - t0
   else
        perturbeEpsion2 = 2.0;
        pertbationMat = sparse(noOfSensors+sDim,noOfSensors+sDim);
        pertbationMat0 = perturbeEpsion2*(ones(noOfSensors,noOfSensors) - noOfSensors*speye(noOfSensors,noOfSensors)); 
        pertbationMat(1:noOfSensors,1:noOfSensors) = pertbationMat0; 
        clear pertbationMat0
    end
    % t0 = cputime;    
    if (pars.objSW == 2)
        c = reshape(pertbationMat,nDim,1);
    else
        % c(nDimLP+1:nDim,1) = reshape(pertbationMat,nDim-nDimLP,1);
        % size(c)
        c = [c(1:nDimLP,:); reshape(pertbationMat,nDim-nDimLP,1)];  
        % size(c)
    end
    % t2 = cputime -t0
    clear pertbationMat
end
% <--- Adding regularization terms to the objective function for anchor free
% cases

%
% fprintf('## cpu time for generating FSDP = %8.2f\n',cputime - startingTime); 
%

% Adding a simplex constraint for the psd variable matrix ---> 
% pars.xAbsBound = 1.0;
% if isfield(pars,'xAbsBound')
%     simplexBound = pars.xAbsBound*(K.s+1);
%     simplexCoefMat = reshape(speye(K.s,K.s),1,K.s*K.s);
%     c = [c(1:K.l,:);sparse(1,1);c(K.l+1:size(A,2),:)];
%     A = [A(:,1:K.l),sparse(size(A,1),1),A(:,K.l+1:size(A,2))];
%     A = [A; [sparse(1,K.l),1,simplexCoefMat]];
%     b = [b; simplexBound];
%     K.l = K.l+1;
% end
% <--- Adding a simplex constraint for the psd variable matrix

% Strengthen the sparsity pattern ---> 
sPatternMat = sparse(K.s,K.s);
sPatternMat(:,K.s-sDim+1:K.s) = 1;%ones(K.s,sDim);
sPatternMat(K.s-sDim+1:K.s,:) = 1;%ones(sDim,K.s);
% if (pars.objSW == 0) || (pars.objSW == 2)
%     sPatternVect = reshape(sPatternMat,1,K.s*K.s);
% else
    sPatternVect = [sparse(1,K.l), reshape(sPatternMat,1,K.s*K.s)];
% end
% <--- Strengthen the sparsity pattern

fprintf('## elapsed time for generating an SDP relaxation problem = %8.2f\n',toc(startingTime));

if  pars.sparseSW == 1
    %%%%%%%%%%
    pars.sDim = sDim;
    %%%%%%%%%%
    
    % Conversion to a sparse LMI SDP ---> 
    [SDP,clique,convMat] = convDual(A,sPatternVect,b,c,K,pars);
    % <--- Conversion to a sparse LMI SDP
    clear A sPatternVect b c
    % Adding upper bounds to \epsilon_{pq}^{+} and \epsilon_{pq}^{-} --->
    % boundSW = 1;
    if (~isempty(pars.noisyFac)) && ((pars.objSW == 1) || (pars.objSW == 3))
        ratio = (2+4*pars.noisyFac)*4*pars.noisyFac;
        distanceVect = reshape(distanceMatrix',1,size(distanceMatrix,1)*size(distanceMatrix,2));
        nzIdxSet = find(distanceVect);
        dimAdd = length(nzIdxSet);
        rowSize = size(SDP.A,1);
        colSize = size(SDP.A,2);
        AMatAdd = sparse(rowSize,dimAdd);
        cAdd = (distanceVect(nzIdxSet) .* distanceVect(nzIdxSet))' * ratio;
        for i=1:dimAdd
            AMatAdd((i-1)*2+1,i) = 1;
            AMatAdd(i*2,i) = 1;
        end
        SDP.A = [SDP.A(:,1:SDP.K.f+SDP.K.l),AMatAdd,SDP.A(:,SDP.K.f+SDP.K.l+1:colSize)];
        SDP.c = [SDP.c(1:SDP.K.f+SDP.K.l,:); cAdd; SDP.c(SDP.K.f+SDP.K.l+1:colSize,:)];
        SDP.K.l = SDP.K.l + dimAdd;
        clear distanceVect AMatAdd
    end
%    clear distanceMatrix
    % <--- Adding upper bounds to \epsilon_{pq}^{+} and \epsilon_{pq}^{-}
    
%    fprintf('## cpu time for generating an SDP relaxation problem = %8.2f\n',cputime - startingTime);

    startingTime = tic; 
    if ~isfield(pars,'SDPsolver') || isempty(pars.SDPsolver) || strcmp(pars.SDPsolver,'sedumi')
        if isnumeric(pars.fid) && (pars.fid == 0)
            fprintf('Starting SeDuMi\n');
        end
        [xbar,ybar,info] = sedumi(SDP.A,SDP.b,SDP.c,SDP.K,pars);
%        info.SDPsolverTime = info.cpusec; 
        if isnumeric(pars.fid) && (pars.fid == 0)
            fprintf('Finished SeDuMi\n');
        end
    elseif strcmp(pars.SDPsolver,'sdpa') % pars.SDPsolver='sdpa;
        OPTION.epsilonStar  = max([pars.eps,1.0e-7]);
        OPTION.epsilonDash  = max([pars.eps,1.0e-7]);
        if (pars.fid == 1)
            OPTION.print = 'display';
        elseif (pars.fid == 0)
            OPTION.print = '';
        elseif isstr(pars.fid)  
            OPTION.print = pars.fid; 
        end
        [xbar,ybar,info] = sedumiwrap(SDP.A,SDP.b,SDP.c,SDP.K,[],OPTION);
%       info.SDPsolverTime = info.sdpaTime; 
    end
    info.SDPsolverTime = toc(startingTime); 
    
    % retriveTime = cputime;
    startingTime = tic; 
    % Retrieving an optimal solution of the original FSDP
    [x] = retrieveFromConvDual(K,ybar,convMat,clique);
    clear convMat clique ybar SDP
    % Extarcting locations of sensors --->
    sdpMatrix = reshape(x(K.l+1:K.l+K.s*K.s,1),K.s,K.s);
    xMatrix = sdpMatrix(noOfSensors+1:K.s,1:noOfSensors);
    xMatrix = [xMatrix, xMatrix0(:,noOfSensors+(1:noOfAnchors))];
    clear K xMatrix0 sdpMatrix
    % <--- Extarcting locations of sensors
%    fprintf('## cpu time for retrieving an optimal solution = %8.2f\n',cputime - retriveTime); 
    fprintf('## elapsed time for retrieving an optimal solution = %8.2f\n',toc(startingTime)); 
    % Retrieving the original ordering for the anchor free case
else % pars.sparseSW == 0
    if (~isempty(pars.noisyFac)) && ((pars.objSW == 1) || (pars.objSW == 3))
        [rowSize,colSize] = size(A); 
        A = [[sparse(rowSize,K.l), A];speye(K.l,K.l),speye(K.l,K.l),sparse(K.l,colSize-K.l)];
        c = [sparse(K.l,1);c];        
        ratio = (2+4*pars.noisyFac)*4*pars.noisyFac;
        distanceVect = reshape(distanceMatrix',1,size(distanceMatrix,1)*size(distanceMatrix,2));
        nzIdxSet = find(distanceVect);
        bAdd = zeros(K.l,1);
        idxSet = 2*[1:length(nzIdxSet)]; 
        bAdd(idxSet,:) = (distanceVect(nzIdxSet) .* distanceVect(nzIdxSet))' * ratio;         
        idxSet = idxSet - 1;
        bAdd(idxSet,:) = (distanceVect(nzIdxSet) .* distanceVect(nzIdxSet))' * ratio;        
        b = [b; bAdd]; 
        K.l = 2*K.l; 
        clear distanceVect bAdd
    end
    
    startingTime = tic;
    if ~isfield(pars,'SDPsolver') || isempty(pars.SDPsolver) || strcmp(pars.SDPsolver,'sedumi')
        if isnumeric(pars.fid) && (pars.fid == 0)
            fprintf('Starting SeDuMi\n');
        end
        [x,y,info] = sedumi(A,b,c,K,pars);
 %       info.SDPsolverTime = info.cpusec; 
        if isnumeric(pars.fid) && (pars.fid == 0)
           fprintf('Finished SeDuMi\n');
        end
    elseif strcmp(pars.SDPsolver,'sdpa') % pars.SDPsolver='sdpa;
        OPTION.epsilonStar  = min([pars.eps,1.0e-7]);
        OPTION.epsilonDash  = min([pars.eps,1.0e-7]);
        if (pars.fid == 1)
            OPTION.print = 'display';
        elseif (pars.fid == 0)
            OPTION.print = '';
        elseif isstr(pars.fid)  
            OPTION.print = pars.fid; 
        end
        [x,y,info] = sedumiwrap(A,b,c,K,[],OPTION);
 %       info.SDPsolverTime = info.sdpaTime; 
    end
    info.SDPsolverTime = toc(startingTime); 

    startingTime = tic;    
    nDim = length(x); 
    x = x(K.l+1:nDim);
    xMatrix = reshape(x,K.s,K.s);
    xMatrix = [xMatrix(K.s-sDim+1:K.s,1:K.s-sDim),xMatrix0(:,noOfSensors+(1:noOfAnchors))];
    clear A b c K xMatrix0
    fprintf('## elapsed time for retrieving an optimal solution = %8.2f\n',toc(startingTime)); 
end

if (noOfAnchors0 == 0) % || (noOfAnchors0 == 1)
    distanceMatrix = [distanceMatrix;distanceMatDeleted];
    distanceMatrix = distanceMatrix + distanceMatrix';
    distanceMatrix = triu(distanceMatrix(perm0,perm0),1);
    xMatrix = xMatrix(:,perm0); 
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [A,b,c,K] = generatePrimalSDP(xMatrix0,noOfAnchors,distanceMatrix,pars)

% 
% This program is based on the Biswas-Ye FSDP, but modified for computational efficiency. 
%
% Input 
% xMatrix0  :   sDim x n matrix of sensors and anchors in the sDim dimensional space, 
%               where n is the number of total sensors and anchors and 
%               anchors are placed in the last m_a columns, 
%               where m_a denotes the number of anchors. 
% noOfAnchors   : the number m_a of anchors; the last m columns of xMatrix0.
% distanceMatrix : the sparse (and noisy) distance matrix between xMatrix0(:,i) and xMatrix0(:,j)

% global ORIGINAL MEMORY
% if MEMORY == 1
%     d = whos('*');
%     mem = 0;
%     for i=1:length(d)
%         mem = mem + d(i).bytes;
%     end
% end

sDim = size(xMatrix0,1); % the dimension of the space
n = size(distanceMatrix,2); % the number of all sensors and anchors

noOfSensors = n-noOfAnchors; % the number of sensors
mDim1 = nnz(distanceMatrix); % the number of distant equations
mDim2 = sDim*(sDim+1)/2; % the number of equalities to specify W = I
mDim = mDim1 + mDim2; 
    % the total number of equalities = 
    % the row size of the constraint matrix A = the number of equality constraitns
if (pars.objSW == 0) || (pars.objSW == 2)
    nDimLP = 0;
else
    nDimLP = 2*mDim1; 
end
    % the number of LP variables; 
    % each nonzero distance yields 2 LP variables
nDimSDP = (noOfSensors+sDim)*(noOfSensors+sDim);
    % the size of the vectorized SDP block in the sedumi format
nDim = nDimLP + nDimSDP;
    % the column size of the constraint matrix A = the number of real
    % variables
K.l = nDimLP;
K.s = noOfSensors+sDim;
%A = sparse(mDim,nDim);
b = sparse(mDim,1);
% cost vector
if (pars.objSW == 1) || (pars.objSW == 3)
    c = sparse((1:nDimLP),1,1,nDim,1,nDimLP);
    %c(1:nDimLP,1) = 1;
else
    c = sparse(nDim,1);
end
% A and b
% if ORIGINAL == 1
%     A = sparse(mDim,nDim);
%     rowPointer = 0;
%     for p = 1:noOfSensors
%         nzIdx = find(distanceMatrix(p,:));
%         for q = nzIdx
%             rowPointer = rowPointer + 1;
%             if (pars.objSW == 1) || (pars.objSW == 3)
%                 uVar = (rowPointer-1)*2 +1;
%                 A(rowPointer,uVar) = 1;
%                 vVar = uVar + 1;
%                 A(rowPointer,vVar) = -1;
%             end
%             if q <= noOfSensors % q is a sensor
%                 Ypp = K.l + (p-1)*K.s +p;
%                 A(rowPointer,Ypp) = 1;
%                 Yqq = K.l + (q-1)*K.s +q;
%                 A(rowPointer,Yqq) = 1;
%                 Ypq = K.l + (p-1)*K.s +q;
%                 A(rowPointer,Ypq) = -1;
%                 Yqp = K.l + (q-1)*K.s +p;
%                 A(rowPointer,Yqp) = -1;
%                 b(rowPointer,1) = distanceMatrix(p,q)*distanceMatrix(p,q);
%             else % q is an anchor q = r
%                 Ypp = K.l + (p-1)*K.s +p;
%                 %fprintf('p = %d, q = %d, rowPointer = %d\n',p,q,rowPointer);
%                 A(rowPointer,Ypp) = 1;
%                 for i=1:sDim
%                     Xip = K.l + (p-1)*K.s + noOfSensors + i;
%                     A(rowPointer,Xip) = -xMatrix0(i,q);
%                     XTip = K.l + (noOfSensors+i-1)*K.s + p;
%                     A(rowPointer,XTip) = -xMatrix0(i,q);
%                 end
%                 b(rowPointer,1) = distanceMatrix(p,q)*distanceMatrix(p,q) - xMatrix0(:,q)'*xMatrix0(:,q);
%                 %fprintf('b(%d) = %18.15e\n',rowPointer,b(rowPointer));
%             end
%         end
%     end
%     for i=1:sDim
%         for j=i:sDim
%             rowPointer = rowPointer + 1;
%             if i == j;
%                 Iii = K.l + (noOfSensors+i-1)*K.s + noOfSensors + i;
%                 A(rowPointer,Iii) = 1;
%                 b(rowPointer,1) = 1;
%             else
%                 Iij = K.l + (noOfSensors+i-1)*K.s + noOfSensors + j;
%                 A(rowPointer,Iij) = 1;
%                 Iji = K.l + (noOfSensors+j-1)*K.s + noOfSensors + i;
%                 A(rowPointer,Iji) = 1;
%                 b(rowPointer,1) = 0;
%             end
%         end
%     end
% elseif ORIGINAL == 0
    if (pars.objSW == 1) || (pars.objSW == 3)
        s = nnz(distanceMatrix(1:noOfSensors,:));
        RowIdx = [(1:s),(1:s)];
        ColIdx = [((1:s)-1)*2+1,((1:s)-1)*2+2];
        Val = [ones(1,s),-ones(1,s)];
    else
        RowIdx = [];
        ColIdx = [];
        Val = [];
    end
    rowPointer = 0;
    for p = 1:noOfSensors
        nzIdx = find(distanceMatrix(p,:));
        if ~isempty(nzIdx)
            sIdx = nzIdx(nzIdx <= noOfSensors);
            aIdx = nzIdx(nzIdx >  noOfSensors);
            if ~isempty(sIdx)
                s = length(sIdx);
                Ypp = K.l + (p-1)*K.s +p;
                Ypp = repmat(Ypp,1,s);
                Yqq = K.l + (sIdx-1)*K.s +sIdx;
                Ypq = K.l + (p-1)*K.s +sIdx;
                Yqp = K.l + (sIdx-1)*K.s +p;

                idx = rowPointer + (1:s);
                RowIdx = [RowIdx,repmat(idx,1,4)];
                ColIdx = [ColIdx,Ypp,Yqq,Ypq,Yqp];
                Val = [Val,ones(1,2*s),-ones(1,2*s)];
                b(idx,1) = distanceMatrix(p,sIdx).*distanceMatrix(p,sIdx);
                rowPointer = rowPointer + s;
            end
            if ~isempty(aIdx)
                Ypp = K.l + (p-1)*K.s +p;
                %disp(aIdx);
                %fprintf('p = %d, rowPointer = %d\n',p,rowPointer);
                idx = repmat(rowPointer+(1:length(aIdx)),2*sDim+1,1);
                idx = idx(:);
                RowIdx = [RowIdx,idx'];
                
                idx = [Ypp,K.l + (p-1)*K.s + noOfSensors + (1:sDim),K.l + (noOfSensors-1+(1:sDim))*K.s + p];
                idx = repmat(idx',1,length(aIdx));
                idx = idx(:);
                ColIdx = [ColIdx,idx'];
                
                tmp = repmat(-xMatrix0(1:sDim,aIdx),2,1);
                tmp = [ones(1,length(aIdx));tmp];
                tmp = tmp(:);
                Val = [Val,tmp'];
                
                %
                % 2008-06-14 Waki
                % I wanted to change the part of making b0 into the part of
                % comment out because the part use for-loop. It consumes
                % much time. 
                %
                % But, the part is a little different from comment out due
                % to numerical errror, so I could not change it. 
                % 
                %{
                tmp = sum(xMatrix0(:,aIdx).*xMatrix0(:,aIdx),1);
                b0(rowPointer+(1:length(aIdx)),1) = distanceMatrix(p,aIdx).*distanceMatrix(p,aIdx) - tmp;
                %}
                
                tmp = zeros(1,length(aIdx));
                for i=1:length(aIdx)
                   tmp(i) =  xMatrix0(:,aIdx(i))'*xMatrix0(:,aIdx(i));
                end
                b(rowPointer+(1:length(aIdx)),1) = distanceMatrix(p,aIdx).*distanceMatrix(p,aIdx) - tmp;
                rowPointer = rowPointer + length(aIdx);
            end
        end
    end
    for i=1:sDim
        for j=i:sDim
            rowPointer = rowPointer + 1;
            if i == j;
                Iii = K.l + (noOfSensors+i-1)*K.s + noOfSensors + i;
                RowIdx = [RowIdx,rowPointer];
                ColIdx = [ColIdx,Iii];
                Val = [Val,1];
                b(rowPointer,1) = 1;
            else
                Iij = K.l + (noOfSensors+i-1)*K.s + noOfSensors + j;
                Iji = K.l + (noOfSensors+j-1)*K.s + noOfSensors + i;
                RowIdx = [RowIdx,rowPointer,rowPointer];
                ColIdx = [ColIdx,Iij,Iji];
                Val = [Val,1,1];
                b(rowPointer,1) = 0;
            end
        end
    end
    %[I,J,V] = find(b-b0);
    %if ~isempty(I)
    %  disp([I';J';V']); 
    %end
    %     size(RowIdx)
    %     size(ColIdx)
    %     size(Val)
    A = sparse(RowIdx,ColIdx,Val,mDim,nDim,length(Val));
% else
%     error('Should set ORIGINAL = 0 or 1.');
% end

% if MEMORY == 1
%     d = whos('*');
%     for i=1:length(d)
%         mem = mem - d(i).bytes;
%     end
%     fprintf('$$ The total Memory in generatePrimalSDP = %4d KB\n',-mem/1024);
% end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [distanceMatrix,distanceMatrixL] = ... 
function [distanceMatrix,noOfNonzeros1d,noOfNonzeros1,maxCliqueSize1] = ... 
    removeRedundantEdges1(sDim,noOfSensors,distanceMatrix,minDegree)
%global ORIGINAL
% distanceMatrixL = distanceMatrix; 
rand('seed',3201); 
[rowSize,colSize] = size(distanceMatrix);
noSensors = rowSize;
%noOfAnchors = colSize - noOfSensors; 
countVector = sparse(1,rowSize); 
for r=noOfSensors+1:colSize
    nzIdx = find(distanceMatrix(:,r)' > 0); 
    if ~isempty(nzIdx)
        for p=nzIdx
            countVector(p) = countVector(p) + 1;
            if (countVector(p) > sDim+1)
                distanceMatrix(p,r) = 0;
%                distanceMatrixL(p,r) = 0;
            end
        end
    end
end

% if ORIGINAL == 1
%     countVector = min([countVector;sDim*ones(1,noOfSensors)],[],1);
%     for q = 1:noSensors
%         idxSensors = find(distanceMatrix(:,q)' > 0);
%         %
%         if ~isempty(idxSensors)
%             %        [temp,idxPermuted] = sort(distanceMatrix(idxSensors,q)');
%             for p=idxSensors
%                 %        for i = idxPermuted
%                 %            p = idxPermuted(i);
%                 %
%                 countVector(p) = countVector(p) + 1;
%                 countVector(q) = countVector(q) + 1;
%                 if (countVector(p) > minDegree) && (countVector(q) > minDegree)
%                     distanceMatrix(p,q) = 0;
%                     %                 if (countVector(p) > 2*minDegree) && (countVector(q) > 2*minDegree)
%                     %                     distanceMatrixL(p,q) = 0;
%                     %                 end
%                 end
%             end
%         end
%     end
%     %full(distanceMatrix)
%     distanceMatrix = sparse(distanceMatrix);
% elseif ORIGINAL == 0

%
    countVector = min([countVector;repmat(sDim,1,noOfSensors)],[],1);
%
    for q = 1:noSensors
        idxSensors = find(distanceMatrix(:,q)' > 0);
        if ~isempty(idxSensors)
            countVector(idxSensors) = countVector(idxSensors)+1;
            s = countVector(q) - minDegree;
            %fprintf('countVector0(%d) = %d, s = %d\n',q,countVector0(q),s);
            idx0 = [];
            if s < 0 && s > -length(idxSensors)
                idx0 = ((1-s):length(idxSensors));
                idx = find(countVector(idxSensors(idx0)) > minDegree);
            elseif s <= -length(idxSensors)
                idx = [];
            else
                idx = find(countVector(idxSensors) > minDegree);
            end
            countVector(q) = countVector(q)+ length(idxSensors);
            if ~isempty(idx)
                if ~isempty(idx0)
                    distanceMatrix(idxSensors(idx0(idx)),q) = 0;
                else
                    distanceMatrix(idxSensors(idx),q) = 0;
                end
            end
            
        end
    end
    %full(distanceMatrix)
    distanceMatrix = sparse(distanceMatrix);
    noOfNonzeros1d = nnz(distanceMatrix);

    sensorDistPatMat = spones(distanceMatrix(:,1:noOfSensors)+distanceMatrix(:,1:noOfSensors)'); 
    sensorDistPatMat = [sensorDistPatMat,ones(noOfSensors,sDim);ones(sDim,noOfSensors),sparse(sDim,sDim)] ...
        + (noOfSensors+1)*speye(noOfSensors+sDim,noOfSensors+sDim);
    permutation = symamd(sensorDistPatMat);
    RMat = chol(sensorDistPatMat(permutation,permutation)); 
    noOfNonzeros1 = nnz(RMat); 
    rVect = sum(spones(RMat),2); 
    maxCliqueSize1 = full(max(rVect'));
    
debugSW = 0;
if debugSW == 1
    sensorDistPatMat = spones(distanceMatrix(:,1:noOfSensors)+distanceMatrix(:,1:noOfSensors)') ...
        + (noOfSensors+1)*speye(noOfSensors,noOfSensors);
    permutation = symrcm(sensorDistPatMat);
    spy(chol(sensorDistPatMat(permutation,permutation)));
    XXXXX
end
%     [I,J,V] = find(distanceMatrix - distanceMatrix);
%     if ~isempty(I)
%         disp(length(I));
%         disp([I';J';V']);
%         disp(full([countVector;countVector0]));
%         fprintf('sDim = %d, noOfSensors = %d\n\n',sDim,noOfSensors);
%         %fprintf('t = %d, t0 = %d\n', t, t0);
%     else
%         
%         fprintf('Same!\n');
%     end
% else
%     error('Should set ORIGINAL = 0 or 1.');
% end
% distanceMatrixL = sparse(distanceMatrixL);
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% end of removeRedundantEdges1 %%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [distanceMatrix,distanceMatrixL] = ... 
function [distanceMatrix,noOfNonzeros2d,noOfNonzeros2,maxCliqueSize2] = ... 
    removeRedundantEdges2(sDim,noOfSensors,distanceMatrix,minDegree)
[rowSize,colSize] = size(distanceMatrix);
% noOfSensors = rowSize;
% noOfAnchors = colSize - noOfSensors; 
% Reordering sensors --->

sensorDistMat = distanceMatrix(:,1:noOfSensors)+distanceMatrix(:,1:noOfSensors)'; 
nodeDegree = full(sum(spones(sensorDistMat),1));
% [temp,permutation] = sort(nodeDegree,'descend');
[temp,permutation] = sort(nodeDegree,'ascend');
[temp,invPermutation] = sort(permutation);
sensorDistMat = sensorDistMat(permutation,permutation); 
distanceMat2 = [sensorDistMat,distanceMatrix(permutation,noOfSensors+1:colSize)];
distanceMat2 = triu(distanceMat2,1);

% <--- Reordering sensors

%global ORIGINAL
% distanceMatrixL = distanceMatrix; 
countVector = sparse(1,rowSize); 
for r=noOfSensors+1:colSize
    nzIdx = find(distanceMat2(:,r)' > 0); 
    if ~isempty(nzIdx)
        for p=nzIdx
            countVector(p) = countVector(p) + 1;
            if (countVector(p) > sDim+1)
                distanceMat2(p,r) = 0;
%                distanceMat2L(p,r) = 0;
            end
        end
    end
end

% if ORIGINAL == 1
%     countVector = min([countVector;sDim*ones(1,noOfSensors)],[],1);
%     for q = 1:noSensors
%         idxSensors = find(distanceMat2(:,q)' > 0);
%         %
%         if ~isempty(idxSensors)
%             %        [temp,idxPermuted] = sort(distanceMat2(idxSensors,q)');
%             for p=idxSensors
%                 %        for i = idxPermuted
%                 %            p = idxPermuted(i);
%                 %
%                 countVector(p) = countVector(p) + 1;
%                 countVector(q) = countVector(q) + 1;
%                 if (countVector(p) > minDegree) && (countVector(q) > minDegree)
%                     distanceMat2(p,q) = 0;
%                     %                 if (countVector(p) > 2*minDegree) && (countVector(q) > 2*minDegree)
%                     %                     distanceMat2L(p,q) = 0;
%                     %                 end
%                 end
%             end
%         end
%     end
%     %full(distanceMat2)
%     distanceMat2 = sparse(distanceMat2);
% elseif ORIGINAL == 0

%
    countVector = min([countVector;repmat(sDim,1,noOfSensors)],[],1);
%

    for q = 1:noOfSensors
        idxSensors = find(distanceMat2(:,q)' > 0);
        if ~isempty(idxSensors)
            countVector(idxSensors) = countVector(idxSensors)+1;
            s = countVector(q) - minDegree;
            %fprintf('countVector0(%d) = %d, s = %d\n',q,countVector0(q),s);
            idx0 = [];
            if s < 0 && s > -length(idxSensors)
                idx0 = ((1-s):length(idxSensors));
                idx = find(countVector(idxSensors(idx0)) > minDegree);
            elseif s <= -length(idxSensors)
                idx = [];
            else
                idx = find(countVector(idxSensors) > minDegree);
            end
            countVector(q) = countVector(q)+ length(idxSensors);
            if ~isempty(idx)
                if ~isempty(idx0)
                    distanceMat2(idxSensors(idx0(idx)),q) = 0;
                else
                    distanceMat2(idxSensors(idx),q) = 0;
                end
            end
            
        end
    end
    %full(distanceMat2)
    distanceMat2 = sparse(distanceMat2);

debugSW = 0;
if debugSW == 1
    sensorDistPatMat = spones(distanceMat2(:,1:noOfSensors)+distanceMat2(:,1:noOfSensors)') ...
        + (noOfSensors+1)*speye(noOfSensors,noOfSensors);
    permutation = symrcm(sensorDistPatMat);
    spy(chol(sensorDistPatMat(permutation,permutation)));
    XXXXX
end
%     [I,J,V] = find(distanceMat2 - distanceMat20);
%     if ~isempty(I)
%         disp(length(I));
%         disp([I';J';V']);
%         disp(full([countVector;countVector0]));
%         fprintf('sDim = %d, noOfSensors = %d\n\n',sDim,noOfSensors);
%         %fprintf('t = %d, t0 = %d\n', t, t0);
%     else
%         
%         fprintf('Same!\n');
%     end
% else
%     error('Should set ORIGINAL = 0 or 1.');
% end
% distanceMat2L = sparse(distanceMat2L);

distanceMat2(:,1:noOfSensors) = distanceMat2(:,1:noOfSensors) + distanceMat2(:,1:noOfSensors)';

debugSW = 0;
if debugSW == 1
    full(distanceMat2)
    countVector
end

distanceMatrix = [distanceMat2(invPermutation,invPermutation),distanceMat2(invPermutation,noOfSensors+1:colSize)];
distanceMatrix = triu(distanceMatrix,1); 
distanceMatrix = sparse(distanceMatrix);

noOfNonzeros2d = nnz(distanceMatrix);

sensorDistPatMat = spones(distanceMatrix(:,1:noOfSensors)+distanceMatrix(:,1:noOfSensors)'); 
sensorDistPatMat = [sensorDistPatMat,ones(noOfSensors,sDim);ones(sDim,noOfSensors),sparse(sDim,sDim)] ...
        + (noOfSensors+1)*speye(noOfSensors+sDim,noOfSensors+sDim);
permutation = symamd(sensorDistPatMat);
RMat = chol(sensorDistPatMat(permutation,permutation));

% spy(RMat+RMat');
% 

noOfNonzeros2 = nnz(RMat);
rVect = sum(spones(RMat),2); 
maxCliqueSize2 = full(max(rVect'));

% XXXXX1127

debugSW = 0;
if debugSW == 1
    sensorDistPatMat = spones(distanceMatrix(:,1:noOfSensors)+distanceMatrix(:,1:noOfSensors)') ...
        + (noOfSensors+1)*speye(noOfSensors,noOfSensors);
    permutation = symamd(sensorDistPatMat);
    spy(chol(sensorDistPatMat(permutation,permutation)));
    XXXXX
end


return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% end of removeRedundantEdges2 %%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [distanceMatrix,distanceMatrixL] = ... 
function [distanceMatrix,noOfNonzeros3d,noOfNonzeros3,maxCliqueSize3] = ... 
    removeRedundantEdges3(sDim,noOfSensors,distanceMatrix,minDegree)
[rowSize,colSize] = size(distanceMatrix);
% noOfSensors = rowSize;
% noOfAnchors = colSize - noOfSensors; 
% Reordering sensors --->

sensorDistMat = distanceMatrix(:,1:noOfSensors)+distanceMatrix(:,1:noOfSensors)'; 
nodeDegree = full(sum(spones(sensorDistMat),1));
% [temp,permutation] = sort(nodeDegree,'descend');
[temp,permutation] = sort(nodeDegree,'descend');
[temp,invPermutation] = sort(permutation);
sensorDistMat = sensorDistMat(permutation,permutation); 
distanceMat2 = [sensorDistMat,distanceMatrix(permutation,noOfSensors+1:colSize)];
distanceMat2 = triu(distanceMat2,1);

% <--- Reordering sensors

%global ORIGINAL
% distanceMatrixL = distanceMatrix; 
countVector = sparse(1,rowSize); 
for r=noOfSensors+1:colSize
    nzIdx = find(distanceMat2(:,r)' > 0); 
    if ~isempty(nzIdx)
        for p=nzIdx
            countVector(p) = countVector(p) + 1;
            if (countVector(p) > sDim+1)
                distanceMat2(p,r) = 0;
%                distanceMat2L(p,r) = 0;
            end
        end
    end
end

% if ORIGINAL == 1
%     countVector = min([countVector;sDim*ones(1,noOfSensors)],[],1);
%     for q = 1:noSensors
%         idxSensors = find(distanceMat2(:,q)' > 0);
%         %
%         if ~isempty(idxSensors)
%             %        [temp,idxPermuted] = sort(distanceMat2(idxSensors,q)');
%             for p=idxSensors
%                 %        for i = idxPermuted
%                 %            p = idxPermuted(i);
%                 %
%                 countVector(p) = countVector(p) + 1;
%                 countVector(q) = countVector(q) + 1;
%                 if (countVector(p) > minDegree) && (countVector(q) > minDegree)
%                     distanceMat2(p,q) = 0;
%                     %                 if (countVector(p) > 2*minDegree) && (countVector(q) > 2*minDegree)
%                     %                     distanceMat2L(p,q) = 0;
%                     %                 end
%                 end
%             end
%         end
%     end
%     %full(distanceMat2)
%     distanceMat2 = sparse(distanceMat2);
% elseif ORIGINAL == 0

%
    countVector = min([countVector;repmat(sDim,1,noOfSensors)],[],1);
%

    for q = 1:noOfSensors
        idxSensors = find(distanceMat2(:,q)' > 0);
        if ~isempty(idxSensors)
            countVector(idxSensors) = countVector(idxSensors)+1;
            s = countVector(q) - minDegree;
            %fprintf('countVector0(%d) = %d, s = %d\n',q,countVector0(q),s);
            idx0 = [];
            if s < 0 && s > -length(idxSensors)
                idx0 = ((1-s):length(idxSensors));
                idx = find(countVector(idxSensors(idx0)) > minDegree);
            elseif s <= -length(idxSensors)
                idx = [];
            else
                idx = find(countVector(idxSensors) > minDegree);
            end
            countVector(q) = countVector(q)+ length(idxSensors);
            if ~isempty(idx)
                if ~isempty(idx0)
                    distanceMat2(idxSensors(idx0(idx)),q) = 0;
                else
                    distanceMat2(idxSensors(idx),q) = 0;
                end
            end
            
        end
    end
    %full(distanceMat2)
    distanceMat2 = sparse(distanceMat2);

debugSW = 0;
if debugSW == 1
    sensorDistPatMat = spones(distanceMat2(:,1:noOfSensors)+distanceMat2(:,1:noOfSensors)') ...
        + (noOfSensors+1)*speye(noOfSensors,noOfSensors);
    permutation = symrcm(sensorDistPatMat);
    spy(chol(sensorDistPatMat(permutation,permutation)));
    XXXXX
end
%     [I,J,V] = find(distanceMat2 - distanceMat20);
%     if ~isempty(I)
%         disp(length(I));
%         disp([I';J';V']);
%         disp(full([countVector;countVector0]));
%         fprintf('sDim = %d, noOfSensors = %d\n\n',sDim,noOfSensors);
%         %fprintf('t = %d, t0 = %d\n', t, t0);
%     else
%         
%         fprintf('Same!\n');
%     end
% else
%     error('Should set ORIGINAL = 0 or 1.');
% end
% distanceMat2L = sparse(distanceMat2L);

distanceMat2(:,1:noOfSensors) = distanceMat2(:,1:noOfSensors) + distanceMat2(:,1:noOfSensors)';

debugSW = 0;
if debugSW == 1
    full(distanceMat2)
    countVector
end

% size(distanceMat2)

distanceMatrix = [distanceMat2(invPermutation,invPermutation),distanceMat2(invPermutation,noOfSensors+1:colSize)];
distanceMatrix = triu(distanceMatrix,1); 
distanceMatrix = sparse(distanceMatrix);
noOfNonzeros3d = nnz(distanceMatrix);

sensorDistPatMat = spones(distanceMatrix(:,1:noOfSensors)+distanceMatrix(:,1:noOfSensors)'); 
sensorDistPatMat = [sensorDistPatMat,ones(noOfSensors,sDim);ones(sDim,noOfSensors),sparse(sDim,sDim)] ...
        + (noOfSensors+1)*speye(noOfSensors+sDim,noOfSensors+sDim);
permutation = symamd(sensorDistPatMat);
RMat = chol(sensorDistPatMat(permutation,permutation));
noOfNonzeros3 = nnz(RMat);
rVect = sum(spones(RMat),2); 
maxCliqueSize3 = full(max(rVect'));

debugSW = 0;
if debugSW == 1
    sensorDistPatMat = spones(distanceMatrix(:,1:noOfSensors)+distanceMatrix(:,1:noOfSensors)') ...
        + (noOfSensors+1)*speye(noOfSensors,noOfSensors);
    permutation = symamd(sensorDistPatMat);
    spy(chol(sensorDistPatMat(permutation,permutation)));
    XXXXX
end


return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% end of removeRedundantEdges3 %%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [distanceMatrix,distanceMatrixL] = ... 
function [distanceMatrix,noOfNonzeros4d,noOfNonzeros4,maxCliqueSize4] = ... 
    removeRedundantEdges4(sDim,noOfSensors,distanceMatrix,minDegree)
[rowSize,colSize] = size(distanceMatrix);
% noOfSensors = rowSize;
% Reordering sensors --->
sensorDistMat = distanceMatrix(:,1:noOfSensors)+distanceMatrix(:,1:noOfSensors)'; 
tempMat = spones(sensorDistMat) + (noOfSensors+1)*speye(noOfSensors,noOfSensors); 
permutation = symrcm(tempMat); 
[tem,invPermutation] = sort(permutation);
sensorDistMat = sensorDistMat(permutation,permutation); 
distanceMat2 = [sensorDistMat,distanceMatrix(permutation,noOfSensors+1:colSize)];
distanceMat2 = triu(distanceMat2,1);
% <--- Reordering sensors
countVector = zeros(1,rowSize);
% Eliminating redundant edges between anchors and sensors ---> 
for r=noOfSensors+1:colSize
    nzIdx = find(distanceMat2(:,r)' > 0); 
    if ~isempty(nzIdx)
        for p=nzIdx
            countVector(p) = countVector(p) + 1;
            if (countVector(p) > sDim+1)
                distanceMat2(p,r) = 0;
            end
        end
    end
end
% <--- Eliminating redundunt edges between anchors and sensors
%
countVector = min([countVector;repmat(sDim,1,noOfSensors)],[],1);
%

DMat = sparse(noOfSensors-1,noOfSensors-1); 
for p=1:noOfSensors-1
    DMat(1:noOfSensors-p,p) = sensorDistMat(p+1:noOfSensors,p);
end

% debugSW = 0;
% if debugSW == 1
%     % 'distanceMatrix'
%     % full(distanceMatrix)
%     % 'sensorDistMat'
%     % full(sensorDistMat)
%     % 'DNat'
%     % full(DMat)
%     % 'countVector'
%     'distanceMat2'
%     full(distanceMat2)
%     countVector
% end

% p = 0;
% nodeDegree = min(countVector);
% while (p <= noOfSensors-2) && (nodeDegree < minDegree)
for p=1:noOfSensors-1    
    nzRowIdx = find(DMat(p,:));
    nzColIdx = nzRowIdx+p;
    countVector(nzRowIdx) = countVector(nzRowIdx) + 1;
    countVector(nzColIdx) = countVector(nzColIdx) + 1;
    edgeDegree = min([countVector(nzRowIdx);countVector(nzColIdx)],[],1);
    idx0 = find(edgeDegree == minDegree+1);
%    debugSW = 0;
%     if debugSW == 1
%         p
%         nzRowIdx
%         nzColIdx
%         countVector
%         edgeDegree
%     end
    if ~isempty(idx0)
%        debugSW = 1;
%         if debugSW == 1
%             full(distanceMat2)
%             idx0
%             nzRowIdx(idx0)
%             nzColIdx(idx0)
%         end
        for q = idx0
            if (countVector(nzRowIdx(q)) > minDegree) && (countVector(nzColIdx(q)) > minDegree)
                distanceMat2(nzRowIdx(q),nzColIdx(q)) = 0;
                countVector(nzRowIdx(q)) = countVector(nzRowIdx(q)) - 1;
                countVector(nzColIdx(q)) = countVector(nzColIdx(q)) - 1;
            end
        end
%        debugSW = 1;
%         if debugSW == 1
%             full(distanceMat2)
%         end
    end
%    debugSW = 1;
%     if (debugSW == 1) && (p==2)
%         XXXXX
%     end
    if min(countVector) >= minDegree
        p = noOfSensors;
    end
end

distanceMat2(:,1:noOfSensors) = distanceMat2(:,1:noOfSensors) + distanceMat2(:,1:noOfSensors)';

debugSW = 0;
if debugSW == 1
    full(distanceMat2)
    countVector
end

% size(distanceMat2)

distanceMatrix = [distanceMat2(invPermutation,invPermutation),distanceMat2(invPermutation,noOfSensors+1:colSize)];
distanceMatrix = triu(distanceMatrix,1); 
distanceMatrix = sparse(distanceMatrix);
noOfNonzeros4d = nnz(distanceMatrix);


sensorDistPatMat = spones(distanceMatrix(:,1:noOfSensors)+distanceMatrix(:,1:noOfSensors)'); 
sensorDistPatMat = [sensorDistPatMat,ones(noOfSensors,sDim);ones(sDim,noOfSensors),sparse(sDim,sDim)] ...
        + (noOfSensors+1)*speye(noOfSensors+sDim,noOfSensors+sDim);
permutation = symamd(sensorDistPatMat);
RMat = chol(sensorDistPatMat(permutation,permutation));
noOfNonzeros4 = nnz(RMat);
rVect = sum(spones(RMat),2); 
maxCliqueSize4 = full(max(rVect'));

debugSW = 0;
if debugSW == 1
    sensorDistPatMat = spones(distanceMatrix(:,1:noOfSensors)+distanceMatrix(:,1:noOfSensors)') ...
        + (noOfSensors+1)*speye(noOfSensors,noOfSensors);
    permutation = symamd(sensorDistPatMat);
    spy(chol(sensorDistPatMat(permutation,permutation)));
    XXXXX
end
% size(distanceMatrix)
% full(distanceMatrix)
% 
% XXXXX
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% end of removeRedundantEdges4 %%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [distanceMatrix,distanceMatrixL] = ... 
function [distanceMatrix,noOfNonzeros5d,noOfNonzeros5,maxCliqueSize5] = ... 
    removeRedundantEdges5(sDim,noOfSensors,distanceMatrix,minDegree)
[rowSize,colSize] = size(distanceMatrix);
% noOfSensors = rowSize;
% noOfAnchors = colSize - noOfSensors; 
% Reordering sensors --->

sensorDistMat = distanceMatrix(:,1:noOfSensors)+distanceMatrix(:,1:noOfSensors)'; 
rand('state',3202); 
nodeDegree = rand(1,noOfSensors);
% [temp,permutation] = sort(nodeDegree,'descend');
[temp,permutation] = sort(nodeDegree);
[temp,invPermutation] = sort(permutation);
sensorDistMat = sensorDistMat(permutation,permutation); 
distanceMat2 = [sensorDistMat,distanceMatrix(permutation,noOfSensors+1:colSize)];
distanceMat2 = triu(distanceMat2,1);

% <--- Reordering sensors

%global ORIGINAL
% distanceMatrixL = distanceMatrix; 
countVector = sparse(1,rowSize); 
for r=noOfSensors+1:colSize
    nzIdx = find(distanceMat2(:,r)' > 0); 
    if ~isempty(nzIdx)
        for p=nzIdx
            countVector(p) = countVector(p) + 1;
            if (countVector(p) > sDim+1)
                distanceMat2(p,r) = 0;
%                distanceMat2L(p,r) = 0;
            end
        end
    end
end

% if ORIGINAL == 1
%     countVector = min([countVector;sDim*ones(1,noOfSensors)],[],1);
%     for q = 1:noSensors
%         idxSensors = find(distanceMat2(:,q)' > 0);
%         %
%         if ~isempty(idxSensors)
%             %        [temp,idxPermuted] = sort(distanceMat2(idxSensors,q)');
%             for p=idxSensors
%                 %        for i = idxPermuted
%                 %            p = idxPermuted(i);
%                 %
%                 countVector(p) = countVector(p) + 1;
%                 countVector(q) = countVector(q) + 1;
%                 if (countVector(p) > minDegree) && (countVector(q) > minDegree)
%                     distanceMat2(p,q) = 0;
%                     %                 if (countVector(p) > 2*minDegree) && (countVector(q) > 2*minDegree)
%                     %                     distanceMat2L(p,q) = 0;
%                     %                 end
%                 end
%             end
%         end
%     end
%     %full(distanceMat2)
%     distanceMat2 = sparse(distanceMat2);
% elseif ORIGINAL == 0

%
    countVector = min([countVector;repmat(sDim,1,noOfSensors)],[],1);
%

    for q = 1:noOfSensors
        idxSensors = find(distanceMat2(:,q)' > 0);
        if ~isempty(idxSensors)
            countVector(idxSensors) = countVector(idxSensors)+1;
            s = countVector(q) - minDegree;
            %fprintf('countVector0(%d) = %d, s = %d\n',q,countVector0(q),s);
            idx0 = [];
            if s < 0 && s > -length(idxSensors)
                idx0 = ((1-s):length(idxSensors));
                idx = find(countVector(idxSensors(idx0)) > minDegree);
            elseif s <= -length(idxSensors)
                idx = [];
            else
                idx = find(countVector(idxSensors) > minDegree);
            end
            countVector(q) = countVector(q)+ length(idxSensors);
            if ~isempty(idx)
                if ~isempty(idx0)
                    distanceMat2(idxSensors(idx0(idx)),q) = 0;
                else
                    distanceMat2(idxSensors(idx),q) = 0;
                end
            end
            
        end
    end
    %full(distanceMat2)
    distanceMat2 = sparse(distanceMat2);

debugSW = 0;
if debugSW == 1
    sensorDistPatMat = spones(distanceMat2(:,1:noOfSensors)+distanceMat2(:,1:noOfSensors)') ...
        + (noOfSensors+1)*speye(noOfSensors,noOfSensors);
    permutation = symrcm(sensorDistPatMat);
    spy(chol(sensorDistPatMat(permutation,permutation)));
    XXXXX
end
%     [I,J,V] = find(distanceMat2 - distanceMat20);
%     if ~isempty(I)
%         disp(length(I));
%         disp([I';J';V']);
%         disp(full([countVector;countVector0]));
%         fprintf('sDim = %d, noOfSensors = %d\n\n',sDim,noOfSensors);
%         %fprintf('t = %d, t0 = %d\n', t, t0);
%     else
%         
%         fprintf('Same!\n');
%     end
% else
%     error('Should set ORIGINAL = 0 or 1.');
% end
% distanceMat2L = sparse(distanceMat2L);

distanceMat2(:,1:noOfSensors) = distanceMat2(:,1:noOfSensors) + distanceMat2(:,1:noOfSensors)';

debugSW = 0;
if debugSW == 1
    full(distanceMat2)
    countVector
end

% size(distanceMat2)

distanceMatrix = [distanceMat2(invPermutation,invPermutation),distanceMat2(invPermutation,noOfSensors+1:colSize)];
distanceMatrix = triu(distanceMatrix,1); 
distanceMatrix = sparse(distanceMatrix);
noOfNonzeros5d = nnz(distanceMatrix);

sensorDistPatMat = spones(distanceMatrix(:,1:noOfSensors)+distanceMatrix(:,1:noOfSensors)'); 
sensorDistPatMat = [sensorDistPatMat,ones(noOfSensors,sDim);ones(sDim,noOfSensors),sparse(sDim,sDim)] ...
        + (noOfSensors+1)*speye(noOfSensors+sDim,noOfSensors+sDim);
permutation = symamd(sensorDistPatMat);
RMat = chol(sensorDistPatMat(permutation,permutation));
noOfNonzeros5 = nnz(RMat);
rVect = sum(spones(RMat),2); 
maxCliqueSize5 = full(max(rVect'));

debugSW = 0;
if debugSW == 1
    sensorDistPatMat = spones(distanceMatrix(:,1:noOfSensors)+distanceMatrix(:,1:noOfSensors)') ...
        + (noOfSensors+1)*speye(noOfSensors,noOfSensors);
    permutation = symamd(sensorDistPatMat);
    spy(chol(sensorDistPatMat(permutation,permutation)));
    XXXXX
end


return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% end of removeRedundantEdges5 %%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [SDP,clique,convMat] = convDual(A,sPatternVect,b,c,K,pars)

% fprintf('## Conversion of a primal SDP into a sparse dual SDP based on the psd completion\n');

%
% Primal SDP ---> Sparse Dual SDP
% 2008/01/30
% Masakazu Kojima
%
% global MEMORY
% if MEMORY == 1
%     d = whos('*');
%     mem = 0;
%     for i=1:length(d)
%         mem = mem + d(i).bytes;
%     end
% end
startingTime = cputime;

SDP = [];
if isfield(K,'q') && ~isempty(K.q)
    fprintf('## Qcone can not be processed\n');
    return 
end
if isfield(K,'r') && ~isempty(K.r)
    fprintf('## Rcone can not be processed\n');
    return
end

if ~isfield(K,'s') || isempty(K.s) 
    fprintf('## no conversion because K.s = []\n'); 
    return
end

% Transpose c if c is a row vector 
if size(c,1) < size(c,2)
    c = c';
end
% Transpose b if b is a row vector 
if size(b,1) < size(b,2)
    b = b';
end

% Check the dimensions of free variables 
if isfield(K,'f') && ~isempty(K.f) && K.f > 0 
    fDim = K.f;
    SDP.K.f = K.f;
else
    fDim = 0; 
    SDP.K.f = 0;
end
% Check the dimensions of LP variables
if isfield(K,'l') && ~isempty(K.l) && K.l > 0 
    ellDim = K.l;
    SDP.K.l = K.l; 
else
    ellDim = 0;
end

%[mDimA,nDimA] = size(A);
mDimA = size(A,1);

% Initialization --->
%    AbarT = []; % To be updated
SDP.A = [];
SDP.b = []; % To be apdated
SDP.c = b; %To be updated, adding 0s
SDP.K.f = mDimA;
SDP.K.s = [];
mDim = 0; % To be updated
nDim = mDimA; % To be updated
colPointer = 0;
% <--- Initialization

% Constraints for free variables
if fDim > 0
    SDP.b = -c(1:fDim,:);
%    AbarT = A(:,1:fDim);
    SDP.A = A(:,1:fDim)'; 
    mDim = fDim;
    colPointer = colPointer + fDim; 
end

% Constraints for LP variables
if ellDim > 0
    SDP.b = [SDP.b;-c(colPointer+(1:ellDim),:)]; 
%    AbarT = [AbarT,A(:,colPointer+1:colPointer+ellDim)];
    SDP.A = [SDP.A;A(:,colPointer+(1:ellDim))'];
    mDim = mDim+ellDim;
    SDP.c = [SDP.c;sparse(ellDim,1)]; 
%    AbarT = [AbarT;[sparse(ellDim,colPointer),-speye(ellDim,ellDim)]]; 
    SDP.A = [SDP.A,[sparse(ellDim,colPointer),-speye(ellDim)]'];
    SDP.K.l = ellDim;
    colPointer = colPointer + ellDim; 
    nDim = nDim + ellDim; 
end 

% Constraints for SDP variables ---> 

% timeForSparsityPattern = 0;
% timeForClique = 0; 
% timeForConversion = 0;

noOfSDPcones = length(K.s);
clique = cell(1,noOfSDPcones);
convMat = cell(1,noOfSDPcones);

% noOfSDPcones

for kk=1:noOfSDPcones
%    startingTime = cputime;
    sDim = K.s(kk);
    if ~isempty(sPatternVect)
        spVect = sPatternVect(colPointer+(1:sDim*sDim));
        if nnz(spVect) == 0
            spVect = [];
        end
    else
        spVect = [];
    end
    Kadd.s = sDim;
    if ~isempty(spVect)
        [sparsityPatternMat] = ...
            genSparsityPatternMat([A(:,colPointer+(1:sDim*sDim));spVect],c(colPointer+(1:sDim*sDim),1),Kadd);
    else
        [sparsityPatternMat] = ...
            genSparsityPatternMat(A(:,colPointer+(1:sDim*sDim)),c(colPointer+(1:sDim*sDim),1),Kadd);
    end
    
%     [rowSize,colSize] = size(sparsityPatternMat);
%     perm = symrcm(sparsityPatternMat(1:rowSize-pars.sDim,1:rowSize-pars.sDim)); 
%     spy(chol(sparsityPatternMat(perm,perm)+rowSize*speye(rowSize-pars.sDim,rowSize-pars.sDim))); 
%     
%     XXXXX

%%%%%%%%%%%%%%%%%%%%%%%% --->
    modifySW = 0;
    if modifySW == 1
%         tempMat = spones(sparsityPatternMat) + (sDim+1)*speye(sDim,sDim);        
%         perm = symamd(tempMat);
%         RMat = chol(tempMat(perm,perm)); 
%         RMat = spones(RMat); 
%         spy(RMat');
%         RVect = (sum(RMat',1));
%         full(RVect(101:200))
%         
%         XXXXX

        perm = symrcm(sparsityPatternMat(1:sDim-2,1:sDim-2));
        
        [temp,invPerm] = sort(perm);
        bandwidth = 15;
        tempMat = sparse(sDim-2,sDim-2);
        for i=1:sDim-2-bandwidth+1
            tempMat(i:i+bandwidth-1,i:i+bandwidth-1) = ones(bandwidth,bandwidth);
        end
        sparsityPatternMat(1:sDim-2,1:sDim-2) = tempMat(invPerm,invPerm);
        perm = symrcm(sparsityPatternMat(1:sDim-2,1:sDim-2));
        spy(sparsityPatternMat);
        %    XXXXXX
    end
%%%%%%%%%%%%%%%%%%%%%%%% --->
    
    %    timeForSparsityPattern = timeForSparsityPattern + cputime - startingTime;
    %    startingTime = cputime;
    %%%%%%%%%%
    % pars.sDim = sDim;
    %%%%%%%%%%
    [oneClique] = cliquesFromSpMatD(sparsityPatternMat,pars);         
    clique{kk} = oneClique;
    %    timeForClique = timeForClique + cputime - startingTime;
    %	startingTime = cputime;
    %     [AbarTAdd,SDP.bAdd,KbarAdd,convMatAdd] = ...
    %         primalToSparseDual(A(:,colPointer+1:colPointer+sDim*sDim),c(colPointer+1:colPointer+sDim*sDim,1),Kadd,clique{kk});
    [AbarAdd,SDP.bAdd,KbarAdd,convMatAdd] = ...
        primalToSparseDual(A(:,colPointer+(1:sDim*sDim)),c(colPointer+(1:sDim*sDim),1),Kadd,clique{kk});
    convMat{kk} = convMatAdd;
    %    timeForConversion = timeForConversion + cputime - startingTime;
    %     nDimAdd = size(AbarTAdd,1) - mDimA;
    %     mDimAdd = size(AbarTAdd,2);
    nDimAdd = size(AbarAdd,2) - mDimA;
    mDimAdd = size(AbarAdd,1);
    % update --->
    %    if isempty(AbarT)
    if isempty(SDP.A)
        %        AbarT = AbarTAdd;
        SDP.A = AbarAdd;
    else
        %         AbarT = [ [AbarT(1:mDimA,:),     AbarTAdd(1:mDimA,:)]; ...
        %             [AbarT(mDimA+1:nDim,:),sparse(nDim-mDimA,mDimAdd)]; ...
        %             [sparse(nDimAdd,mDim), AbarTAdd(mDimA+1:mDimA+nDimAdd,:)] ];
        SDP.A = [ [SDP.A(:,1:mDimA); AbarAdd(:,1:mDimA)], ...
            [SDP.A(:,mDimA+1:nDim); sparse(mDimAdd,nDim-mDimA)], ...
            [sparse(mDim,nDimAdd); AbarAdd(:,mDimA+(1:nDimAdd))] ];
    end
    SDP.b = [SDP.b; SDP.bAdd];
    SDP.c = [SDP.c; sparse(nDimAdd,1)];
    SDP.K.s = [SDP.K.s, KbarAdd.s];
    mDim = mDim + mDimAdd;
    nDim = nDim + nDimAdd;
    % <--- update
    colPointer = colPointer + sDim*sDim;
end
% <--- Constraints for SDP variables

% size(SDP.A)

% fprintf('## cpu time for computing a sparsity pattern = %6.2e\n',timeForSparsityPattern);
% fprintf('## cpu time for computing cliques = %6.2e\n',timeForClique);

%
% fprintf('## cpu time for conversion = %8.2f\n',cputime - startingTime);
%

debugSW = 0;
if debugSW == 1
    %    Abar = AbarT';
    for kk=1:noOfSDPcones
        fprintf('clique{%d}.Set\n',kk);
        full(clique{kk}.Set)
    end
    fprintf('SDP.K.s = \n');
    SDP.K
    SDP.K.s
    noOfSDPcones = length(SDP.K.s);
    %    [rowSize,colSize] = size(Abar)
    %[rowSize,colSize] = size(SDP.A)
    rowSize = size(SDP.A,1);
    fprintf('Cbar = \n');
    pointer = SDP.K.f;
    for k=1:noOfSDPcones
        full(reshape(SDP.c(pointer+1:pointer+SDP.K.s(k)*SDP.K.s(k),1),SDP.K.s(k),SDP.K.s(k)))
        pointer = pointer + SDP.K.s(k)*SDP.K.s(k);
    end
    for i=1:rowSize
        fprintf('SDP.A{%d} = \n',i);
        pointer = SDP.K.f;
        for k=1:noOfSDPcones
            full(reshape(SDP.A(i,pointer+1:pointer+SDP.K.s(k)*SDP.K.s(k)),SDP.K.s(k),SDP.K.s(k)))
            pointer = pointer + SDP.K.s(k)*SDP.K.s(k);
        end
    end
    XXXXX
end
% SDP.A = AbarT';

% if MEMORY == 1
%     d = whos('*');
%     for i=1:length(d)
%         mem = mem - d(i).bytes;
%     end
%     fprintf('$$ The total Memory in convDual = %4d KB\n',-mem/1024);
% end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [clique] = cliquesFromSpMatD(sparsityPatternMat,pars)
%
% 2008-06-13 Waki
% Caution!
% I have changed the structures of clique.
% 
% Before:
% the structures of clique are 
%     NoC, maxC, minC, Set and idxMatrix.
%
% After:
% the structures of clique are 
%     NoC, maxC, minC, Elem, NoElem and idxMatrix.
%
% Elem is a row vector, which contains all cliques.
% NoElem is NOC-dimensional row vector, which has each size of cliques
%
% For example, we consider three cliques
% c1 ={1,2}, c2 = {1,3,4}, c3 ={5}.
%
% Then 
% clique.Elem =[1,2,1,3,4,5]
% clique.NoElem = [2,3,1]
%
% If you want to access the i-th clqiue, execute the following command:
%
% idx = sum(clique.NoElem(1:i-1));
% clique.Elem(idx+(1:clique.NoElem(i)));
%
% In this example, if you want to access the second clique, you should execute
%
% idx = clique.NoElem(1);
% clique.Elem(idx+(1:clique.NoElem(2)));
%

% global ORIGINAL MEMORY
% 
% if MEMORY == 1
%     d = whos('*');
%     mem = 0;
%     for i=1:length(d)
%         mem = mem + d(i).bytes;
%     end
% end

%spones(sparsityPatternMat);
nDim = size(sparsityPatternMat,1);
sparsityPatternMat = spones(sparsityPatternMat) + (2*nDim+1)*speye(nDim);
%%%%%
orderingSW = 0;
if orderingSW == 0%
% pars.edgeRemoveSW == 0;
    %% minimum degree ordering
%    orderingSW = 0;
    I = symamd(sparsityPatternMat);
elseif orderingSW == 1
    %% sparse reverse Cuthill-McKee ordering
%    orderingSW = 1;
    I = symrcm(sparsityPatternMat);   
elseif orderingSW == 3
% pars.edgeRemoveSW = 3;
    %% sparse reverse Cuthill-McKee ordering only for the sensor network
    %% problem
%    orderingSW = 1;
    I = symrcm(sparsityPatternMat(1:nDim-pars.sDim,1:nDim-pars.sDim));
    I = [I,(nDim-pars.sDim+1):nDim];     
%    I = symamd(sparsityPatternMat);
% if orderingSW == 0
%     
%     %% minimum degree ordering
%     I = symamd(sparsityPatternMat);
% elseif orderingSW == 1
%     %% sparse reverse Cuthill-McKee ordering
%     I = symrcm(sparsityPatternMat);
end
%% cholesky decomposition
[R,p] = chol(sparsityPatternMat(I,I));
if (p > 0)
    error('Correlative sparsity matrix is not positive definite.');
end

debug = 0;
if debug == 1
    RR = R+R';
    spy(RR);
    XXXXX
end
%%
%% Step3
%% Finding the maxmal clieques
%%
%% put 1 for nonzero element of R
% if ORIGINAL == 1
%     Cliques = spones(R);
%     [value,orig_idx] = sort(I);
%     remainIdx = 1;
%     for i=2:nDim
%         checkSet = Cliques(i,i:nDim);
%         one = find(checkSet);
%         noOfone = length(one);
%         cliqueResult = Cliques(1:i-1,i:nDim)*checkSet';
%         yesno = find(cliqueResult == noOfone);
%         %%
%         %% Remove the set included into other set.
%         %%
%         if ~any(yesno)
%             remainIdx = [remainIdx;i];
%         end
%     end
%     clique.Set = Cliques(remainIdx,orig_idx);
%     %%
%     %% Clique Information
%     %%
%     clique.NoC  = size(clique.Set,1);
%     sumClique = full(sum(clique.Set,2));
%     clique.maxC = full(max(sumClique));
%     clique.minC = full(min(sumClique));
%     %% Information on indexing variables
%     % Cliques = spones(R+R');
%     % Cliques = Cliques(orig_idx,orig_idx);
%     Cliques = sparse(nDim,nDim);
%     for i=1:clique.NoC
%         idx = find(clique.Set(i,:));
%         sDimE = length(idx);
%         Cliques(idx,idx) = ones(sDimE,sDimE);
%     end
%     Cliques = triu(Cliques);
%     % % spy(Cliques);
%     pointer = 0;
%     for i=1:nDim
%         nnzRowIdx = find(Cliques(i,:));
%         noNnz = length(nnzRowIdx);
%         Cliques(i,nnzRowIdx) = [pointer+1:pointer+noNnz];
%         pointer = pointer + noNnz;
%     end
%     % full(Cliques)
%     clique.idxMatrix = Cliques;
% elseif ORIGINAL == 0
    Cliques = spones(R);
    [value,orig_idx] = sort(I);
    remainIdx = 1;
    for i=2:nDim
        idx = i:nDim;
        one = find(Cliques(i,idx));
        noOfone = length(one);
        %
        % 2008-06-08 Waki
        % replace multiplication by sum.
        %
        cliqueResult = sum(Cliques(remainIdx,idx(one)),2);
        if isempty(find(full(cliqueResult) == noOfone,1))
            remainIdx = [remainIdx;i];
        end
    end
    cSet = Cliques(remainIdx,orig_idx);
    %%
    %% Clique Information
    %%
    clique.NoC  = length(remainIdx);
    [I,J] = find(cSet');
    clique.Elem = I;
    clique.NoElem = full(sum(cSet,2));
    clique.NoElem = clique.NoElem';
    clique.maxC = full(max(clique.NoElem));
    clique.minC = full(min(clique.NoElem));
    Cliques = sparse(nDim,nDim);
    idx = 0;
    for i=1:clique.NoC
        s = clique.NoElem(i);
        tmp = clique.Elem(idx+(1:s));
        idx = idx + s;
        %
        % 2008-06-08 Waki
        % replace ones by 1
        %
        Cliques(tmp,tmp) = 1;
    end
    Cliques = tril(Cliques);
    
    %
    % 2008-06-10 Waki
    % Replace the part of renumbering Cliques by a faster version.
    %
    [I,J,V] = find(Cliques);
    s = length(V);
    clique.idxMatrix = sparse(J,I,(1:s),nDim,nDim,s);
% else
%     error('Should set ORIGINAL = 0 or 1.');
% end
% 
% if MEMORY == 1
%     d = whos('*');
%     for i=1:length(d)
%         mem = mem - d(i).bytes;
%     end
%     fprintf('$$ The total Memory in cliquesFromSpMatD = %4d KB\n',-mem/1024);
% end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Abar,bbar,Kbar,convMat] = primalToSparseDual(A,c,K,clique)

% global ORIGINAL MEMORY
% if MEMORY == 1
%     d = whos('*');
%     mem = 0;
%     for i=1:length(d)
%         mem = mem + d(i).bytes;
%     end
% end

%
% 20080-06-07 Waki 
% Replace this part by a faster version.
%
% if ORIGINAL == 1
%     % Equality constraint matrix in dual format --->
%     %AbarT = [];
%     Abar = [];
%     bbar = [];
%     sDim = K.s;
%     for row = 1:sDim
%         for col = find(clique.idxMatrix(row,:));
%             posInA = (row-1)*sDim+col;
%             if col == row
%                 bbar = [bbar;-c(posInA,:)];
%                 %            AbarT = [AbarT,A(:,posInA)];
%                 Abar = [Abar; A(:,posInA)'];
%             else
%                 bbar = [bbar; -2*c(posInA,:)];
%                 %            AbarT = [AbarT,2*A(:,posInA)];
%                 Abar = [Abar; 2*A(:,posInA)'];
%             end
%         end
%     end
% elseif ORIGINAL == 0
    if size(c,1) > 1 && size(c,2) > 1
       error('c must be a vector, not a matrix.'); 
    end
    sDim = K.s;
    [J,I] = find(clique.idxMatrix');
    posInA = (I-1)*sDim + J;
    bbar = -2*c(posInA);
    Abar = 2*A(:,posInA);
    eqIdx = find(J == I);
    bbar(eqIdx) = -c(posInA(eqIdx));
    Abar(:,eqIdx) = A(:,posInA(eqIdx));
    Abar = Abar';
% else
%     error('Should set ORIGINAL = 0 or 1.');
% end
% [rowSizeAbarT,colSizeAbarT] = size(AbarT);

%
% 2008-06-07 Waki
% replace this part by a faster version.
%
% The part of "AbarSDPT = [AbarSDPT; -convMat{i}];" is faster than the
% version where we allocate AbarSDPT.
%
%
% if ORIGINAL == 1
%     [rowSizeAbarT,colSizeAbarT] = size(Abar');
%     % <--- Equality constraint matrix in dual format
%     Kbar.s = [];
%     sdpDim = colSizeAbarT;
% 
%     % Sparse SDP constraint matrix in dual format --->
%     for i=1:clique.NoC
%         sDimE = nnz(clique.Set(i,:));
%         Kbar.s = [Kbar.s,sDimE];
%         convMat{i} = sparse(sDimE*sDimE,sdpDim);
%     end
%     for i=1:clique.maxC
%         psdConstMat{i} = [];
%     end
% 
%     AbarSDPT = [];
%     for i=1:clique.NoC
%         idxSet = find(clique.Set(i,:));
%         sDimE = length(idxSet);
%         if isempty(psdConstMat{sDimE})
%             [psdMat] = genPsdConstMat(sDimE);
%             psdConstMat{sDimE} = psdMat;
%         else
%             psdMat = psdConstMat{sDimE};
%         end
%         tempVector = reshape(clique.idxMatrix(idxSet,idxSet)',1,sDimE*sDimE);
%         tempIdx = find(tempVector);
%         idxInAbarSDPT = tempVector(tempIdx);
%         convMat{i}(:,idxInAbarSDPT) = psdMat;
%         AbarSDPT = [AbarSDPT; -convMat{i}];
%     end
%     % <--- Sparse SDP constraint matrix in dual format
%     % AbarT = [AbarT; AbarSDPT];
%     Abar = [Abar, AbarSDPT'];
% elseif ORIGINAL == 0
    [rowSizeAbarT,sdpDim] = size(Abar'); 
    % <--- Equality constraint matrix in dual format 

    % Sparse SDP constraint matrix in dual format --->
    psdConstMat = cell(1,clique.maxC);
    for i=1:clique.maxC
        psdConstMat{i} = [];
    end
    convMat = cell(1,clique.NoC);
    Kbar.s = clique.NoElem;
    idx = 0;
    for i=1:clique.NoC
        sDimE = clique.NoElem(i);
        tmpIdx = clique.Elem(idx+(1:sDimE));
        idx = idx + sDimE;
        tempVector = clique.idxMatrix(tmpIdx,tmpIdx)';
        tempVector = tempVector(:);
        tempVector = tempVector';
        idxInAbarSDPT = tempVector(tempVector ~= 0);
        if isempty(psdConstMat{sDimE})
            [RowIdx, ColIdx] = genPsdConstMat2(sDimE);
            convMat{i} = sparse(idxInAbarSDPT(ColIdx),RowIdx,1,sdpDim,sDimE*sDimE,sDimE*sDimE);
        else
            psdMat = psdConstMat{sDimE};
            convMat{i}(idxInAbarSDPT,:) = psdMat';
        end
    end
    AbarSDPT = [convMat{1:clique.NoC}];
    % <--- Sparse SDP constraint matrix in dual format
    % AbarT = [AbarT; AbarSDPT];
    Abar = [Abar, -AbarSDPT];
% else
%     error('Should set ORIGINAL = 0 or 1.');
% end
% 
% if MEMORY == 1
%     d = whos('*');
%     for i=1:length(d)
%         mem = mem - d(i).bytes;
%     end
%     fprintf('$$ The total Memory in primalToSparseDual = %4d KB\n',-mem/1024);
% end

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [psdMat] = genPsdConstMat(sDimE)
%global ORIGINAL
%
% 2008-06-07 Waki
% Replace this part by a faster version.
%
% if ORIGINAL == 1
%     psdMat = [];
%     % size(psdMat)
%     % sDimE
%     for row = 1:sDimE
%         for col = row:sDimE
%             oneColumn = sparse(sDimE*sDimE,1);
%             oneColumn((row-1)*sDimE+col,1) = 1;
%             if row < col
%                 oneColumn((col-1)*sDimE+row,1) = 1;
%             end
%             %         size(psdMat)
%             %         size(oneColumn)
%             psdMat = [psdMat,oneColumn];
%         end
%     end
% elseif ORIGINAL == 0
    RowIdx = zeros(sDimE*sDimE,1);
    ColIdx = zeros(sDimE*sDimE,1);
    idx = 1;
    k = 1;
    for row = 1:sDimE
        RowIdx(idx) = (row-1)*sDimE + row;
        ColIdx(idx) = k;
        idx = idx + 1;
        k = k + 1;
        for col = (row+1):sDimE
            RowIdx(idx) = (row-1)*sDimE + col;
            ColIdx(idx) = k;
            idx = idx + 1;
            RowIdx(idx) = (col-1)*sDimE + row;
            ColIdx(idx) = k;
            idx = idx + 1;
            k = k + 1;
        end
    end
    psdMat = sparse(RowIdx,ColIdx,1,sDimE*sDimE,sDimE*(sDimE+1)/2,sDimE*sDimE);
    clear RowIdx ColIdx
% else
%     error('Should set ORIGINAL = 0 or 1.');
% end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 2008-06-08 Waki 
% Revise genPsdConstMat. Change outputed variable.
%
%
function [RowIdx, ColIdx] = genPsdConstMat2(sDimE)

[I,J] = find(ones(sDimE));
RowIdx = (J-1)*sDimE + I;
tmpMat = zeros(sDimE);
idx = 0;
for i=1:sDimE
    tmpMat(i,i:sDimE) = idx + (i:sDimE);
    idx = idx + sDimE -i;
end
ColIdx = triu(tmpMat)+triu(tmpMat,1)';
ColIdx = ColIdx(:);
clear tmpMat
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sparsityPatternMat] = genSparsityPatternMat(A,c,K)
% global MEMORY
% 
% if MEMORY == 1
%     d = whos('*');
%     mem = 0;
%     for i=1:length(d)
%         mem = mem + d(i).bytes;
%     end
% end
sparsityPatternMat = [];
% if isfield(K,'q') && ~isempty(K.q) && (K.q > 0) 
%     fprintf('## Qcone (K.q > 0) can not be processed.\n');
%     return 
% end
% if isfield(K,'r') && ~isempty(K.r) &&  (K.r > 0) 
%     fprintf('## Rcone (K.r > 0) can not be processed.\n');
%     return
% end
if ~isfield(K,'s') || isempty(K.s)  
    fprintf('## K.s is not assigned\n');
    return
elseif length(K.s) > 1
    fprintf('## length(K.s) = 1 is assumed in this implementation\n');
    return
end

fDim = 0;
if isfield(K,'f') && ~isempty(K.f) && K.f > 0 
    fDim = K.f;  
end

ellDim = 0;
if isfield(K,'l') && ~isempty(K.l) && K.l > 0 
    ellDim = K.l;
end

qDim = 0;
if isfield(K,'q') && ~isempty(K.q)
    qDim = sum(K.q); 
end

rDim = 0;
if isfield(K,'r') && ~isempty(K.r)
    rDim = sum(K.r); 
end

nonSDPdim = fDim+ellDim+qDim+rDim; 

%mDim = size(A,1); 
% Dimensions of SDP variables
sDim = 0; 
if isfield(K,'s') && ~isempty(K.s)
    sDim = sum(K.s);
end

if size(c,2) > size(c,1)
    c = c';
end
% Construction of sparsityPatternMat
sparsityPatternMat = speye(sDim);

if isfield(K,'s') && ~isempty(K.s)
    pointer = 0;
    for p=1:length(K.s)
        kDim = K.s(p); 
        tempVec = abs(c(nonSDPdim+pointer+(1:kDim*kDim),1))'; 
        tempVec = tempVec + sum(abs(A(:,nonSDPdim+pointer+(1:kDim*kDim))),1);
        tempMat = reshape(tempVec,kDim,kDim); 
        sparsityPatternMat(pointer+(1:kDim),pointer+(1:kDim)) = sparsityPatternMat(pointer+(1:kDim),pointer+(1:kDim)) + tempMat; 
        pointer = pointer + kDim; 
    end
end
sparsityPatternMat = spones(sparsityPatternMat);

% spy(sparsityPatternMat); 

% if MEMORY == 1
%     d = whos('*');
%     for i=1:length(d)
%         mem = mem - d(i).bytes;
%     end
%     fprintf('$$ The total Memory in genSparsityPatternMat = %4d KB\n',-mem/1024);
% end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xPrimal] = retrieveFromConvDual(K0,ybar,convMat,clique)
% global ORIGINAL MEMORY
% if MEMORY == 1
%     d = whos('*');
%     mem = 0;
%     for i=1:length(d)
%         mem = mem + d(i).bytes;
%     end
% end

if isfield(K0,'f') && ~isempty(K0.f) && K0.f > 0
    fDim = K0.f;
else
    fDim = 0; 
end
if isfield(K0,'l') && ~isempty(K0.l) && K0.l > 0 
    ellDim = K0.l;
else
    ellDim = 0; 
end

if fDim > 0
    xPrimalFree = ybar(1:fDim,1); 
else
    xPrimalFree =[];
end

if ellDim > 0
    xPrimalLP = ybar(fDim+1:fDim+ellDim,1); 
else
    xPrimalLP =[];
end

noOfSDPcones = length(K0.s);
ybarPointer = fDim+ellDim;
xPrimalSDP = [];
for kk=1:noOfSDPcones
    sDim = K0.s(kk);
    xAddMat = sparse(sDim,sDim);
    % size(xPrimalSDP)
%     if ORIGINAL == 1
%         kDim = size(convMat{kk}{1},2);
%     elseif ORIGINAL == 0
        kDim = size(convMat{kk}{1},1);
%     else
%         error('Should set ORIGINAL = 0 or 1.');
%     end
    idx = 0;
    for p=1:clique{kk}.NoC
%         if ORIGINAL == 1
%             pVect = convMat{kk}{p} * ybar(ybarPointer+(1:kDim));
%             cliqueIdx = find(clique{kk}.Set(p,:));
%             sDimE = length(cliqueIdx);
%             pMat = reshape(pVect,sDimE,sDimE);
%             xAddMat(cliqueIdx,cliqueIdx) = pMat;
%         elseif ORIGINAL == 0
            %
            % 2008-06-10 Waki
            % Replace cliqueIdx by clique.elem{i}. This variable is defined in
            % cliquesFromSpMatD.
            %
            pVect = convMat{kk}{p}' * ybar(ybarPointer+(1:kDim));
            sDimE = clique{kk}.NoElem(p);
            pMat = reshape(pVect,sDimE,sDimE);
            tmpIdx = clique{kk}.Elem(idx+(1:sDimE));
            xAddMat(tmpIdx,tmpIdx) = pMat;
            idx = idx +sDimE;
%         else
%             error('Should set ORIGINAL = 0 or 1.');
%         end
    end
    ybarPointer = ybarPointer+kDim; 
    xPrimalSDP = [xPrimalSDP; reshape(xAddMat,sDim*sDim,1)]; 
end

xPrimal = [xPrimalFree; xPrimalLP; xPrimalSDP];

% if MEMORY == 1
%     d = whos('*');
%     for i=1:length(d)
%         mem = mem - d(i).bytes;
%     end
%     fprintf('$$ The total Memory in retrieveFromConvDual = %4d KB\n',-mem/1024);
% end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [clique,anchorLocation] = findClique(distanceMatrix,cliqueSize);

clique = [];
anchorLocation = [];
[noOfSensors,colSize] = size(distanceMatrix); 
if cliqueSize < 3
    fprintf('# cliqueSize = %d < 3 = %d\n',cliqueSize);
    return
elseif 4 < cliqueSize
    fprintf('# cliqueSize = %d > 4 = %d\n',cliqueSize);
    return
elseif noOfSensors < colSize
    fprintf('# noOfSensors = %d < colSize = %d\n',noOfSensors,colSize);
    return
elseif noOfSensors > colSize
    fprintf('# noOfSensors = %d > colSize = %d\n',noOfSensors,colSize);
    return
elseif noOfSensors < 3
    fprintf('# noOfSensors = %d < cliqueSize = %d\n',noOfSensors,cliqueSize);
    return    
end

incidenceMatrix = spones(distanceMatrix+distanceMatrix')+speye(noOfSensors,noOfSensors); 

if cliqueSize == 3 
    [clique] = find3Clique(incidenceMatrix);
%    clique
    i1 = clique(1); x1 = zeros(2,1); 
    i2 = clique(2); x2 = zeros(2,1); 
    i3 = clique(3); x3 = zeros(2,1); 
    x2(1,1) = distanceMatrix(i1,i2);
    x3(1,1) = (distanceMatrix(i1,i3)^2-distanceMatrix(i2,i3)^2+x2(1,1)^2)/(2*x2(1,1)); 
    a = distanceMatrix(i1,i3)^2 - x3(1,1)^2; 
    if a < 0
%        fprintf('## Some error computing 3 points\n');
        anchorLocation = [];
        return
    end
    x3(2,1) = sqrt(a);
    anchorLocation = [x1,x2,x3];
elseif cliqueSize == 4 
    [clique] = find4Clique(incidenceMatrix); 
    i1 = clique(1); x1 = zeros(3,1);
    i2 = clique(2); x2 = zeros(3,1);
    i3 = clique(3); x3 = zeros(3,1);
    i4 = clique(4); x4 = zeros(3,1);
    x2(1,1) = distanceMatrix(i1,i2);
    x3(1,1) = (distanceMatrix(i1,i3)^2-distanceMatrix(i2,i3)^2+x2(1,1)^2)/(2*x2(1,1)); 
    a = distanceMatrix(i1,i3)^2 - x3(1,1)^2; 
    if a < 0
%        fprintf('# Some error computing 3 points\n');
%         x2(1,1)
%         x3(1,1)
%         x3(2,1) = 0.0;
        anchorLocation = [];
%         clique
%         full(distanceMatrix(clique,clique))
%        XXXX
        return
    else
        x3(2,1) = sqrt(a);        
    end
    bVect = [distanceMatrix(i2,i4)^2-distanceMatrix(i1,i4)^2-x2'*x2;distanceMatrix(i3,i4)^2-distanceMatrix(i1,i4)^2-x3'*x3];
    AMat = -2*[x2(1,1), 0; x3(1,1), x3(2,1)];
    y = AMat \ bVect;
    x4(1,1) = y(1,1);
    x4(2,1) = y(2,1); 
    b = distanceMatrix(i1,i4)^2 -x4(1,1)^2 -x4(2,1)^2; 
    if b < 0
%        fprintf('# Some error computing 4 points\n');
        x4(3,1) = 0.0;
        anchorLocation = [];
        return
    end
    x4(3,1) = sqrt(b); 
	anchorLocation = [x1,x2,x3,x4];
else
    return
end

debugSW = 0;
if debugSW == 1
    format long
    mDim = length(clique); 
    clique
    full(distanceMatrix(clique,clique))
    dMat = zeros(mDim,mDim);
    for i=1:mDim
        p=clique(i);
        for j=i+1:mDim
            q=clique(j);
            dMat(i,j) = norm(anchorLocation(:,i)'-anchorLocation(:,j)');
        end
    end
    full(dMat)
    format short
    XXXXX
end

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cliqueOut] = checkclique(incidenceMatrix,cliqueIn) 

% incidenceMatrix --- nodes to nodes incident matrix
% 	0, 1 matrix
%   = 1 if i = j
%   = 1 if the node i is incident to tne node j
%   = 0 otherwise

cliqueOut = [];
if isempty(cliqueIn) 
    return
elseif length(unique(cliqueIn)) < length(cliqueIn) 
    return
end
kDim = length(cliqueIn); 
checkVector = reshape(incidenceMatrix(cliqueIn,cliqueIn),1,kDim*kDim); 
if isempty(find(checkVector == 0)) 
    cliqueOut = cliqueIn;
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [clique] = find3Clique(incidenceMatrix)

clique = [];
[noOfSensors,colSize] = size(incidenceMatrix); 
if noOfSensors < colSize
    fprintf('noOfSensors = %d < colSize = %d\n',noOfSensors,colSize);
    return
elseif noOfSensors > colSize
    fprintf('noOfSensors = %d > colSize = %d\n',noOfSensors,colSize);
    return
elseif noOfSensors < 3
    fprintf('noOfSensors = %d < 3\n',noOfSensors);
    return    
end

%%%%%
% for debug
% rand('state',3201); 
% incidenceMatrix = spones(incidenceMatrix+incidenceMatrix'+speye(noOfSensors,noOfSensors));
% [temp,perm] = sort(rand(1,noOfSensors));
% incidenceMatrix = incidenceMatrix(perm,perm);
% full(incidenceMatrix)
%%%%%

p=0;
controlSW = 0;
while controlSW == 0
    p = p+1;
    pthRowNzIdx = find(incidenceMatrix(p,p+1:noOfSensors));
    if ~isempty(pthRowNzIdx)
        pthRowNzIdx = pthRowNzIdx + p;
        pDim = length(pthRowNzIdx); 
        i = 0;
        while (i < pDim) && (controlSW == 0)
            i = i+1;
            q = pthRowNzIdx(i); 
            if q < noOfSensors
                qthRowNzIdx = find(incidenceMatrix(q,q+1:noOfSensors));
                if ~isempty(qthRowNzIdx)
                    qthRowNzIdx = qthRowNzIdx + q;
                    qDim = length(qthRowNzIdx); 
                    j = 0;
                    while (j < qDim) && (controlSW == 0)
                        j = j+1;
                        r = qthRowNzIdx(j);
                        [clique] = checkclique(incidenceMatrix,[p,q,r]);
                        if ~isempty(clique)
                            controlSW = 1;
                        end
                    end
                end
            end
        end
    end
    if (p == noOfSensors-2) && (controlSW == 0)
        controlSW = -1;
    end
end

% clique

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [clique] = find4Clique(incidenceMatrix)

clique = [];
[noOfSensors,colSize] = size(incidenceMatrix); 
if noOfSensors < colSize
    fprintf('noOfSensors = %d < colSize = %d\n',noOfSensors,colSize);
    return
elseif noOfSensors > colSize
    fprintf('noOfSensors = %d > colSize = %d\n',noOfSensors,colSize);
    return
elseif noOfSensors < 4
    fprintf('noOfSensors = %d < 4\n',noOfSensors);
    return    
end

%%%%%
% for debug
% rand('state',3204); 
% incidenceMatrix = spones(incidenceMatrix+incidenceMatrix'+speye(noOfSensors,noOfSensors));
% [temp,perm] = sort(rand(1,noOfSensors));
% incidenceMatrix = incidenceMatrix(perm,perm);
% full(incidenceMatrix)
%%%%%

% [clique] = checkclique(incidenceMatrix,[1,3,8,9])

p=0;
controlSW = 0;
while controlSW == 0
    p = p+1;
    pthRowNzIdx = find(incidenceMatrix(p,p+1:noOfSensors));
    if ~isempty(pthRowNzIdx)
        pthRowNzIdx = pthRowNzIdx + p;
        pDim = length(pthRowNzIdx); 
        i = 0;
        while (i < pDim) && (controlSW == 0)
            i = i+1;
            q = pthRowNzIdx(i); 
            if q < noOfSensors-1
                qthRowNzIdx = find(incidenceMatrix(q,q+1:noOfSensors));
                if ~isempty(qthRowNzIdx)
                    qthRowNzIdx = qthRowNzIdx + q;
                    qDim = length(qthRowNzIdx); 
                    j = 0;
                    while (j < qDim) && (controlSW == 0)
                        j = j+1;
                        r = qthRowNzIdx(j);
                        if r < noOfSensors 
                            rthRowNzIdx = find(incidenceMatrix(r,r+1:noOfSensors));
                            if ~isempty(rthRowNzIdx)
                                rthRowNzIdx = rthRowNzIdx + r;
                                rDim = length(rthRowNzIdx); 
                                k = 0; 
                                while (k < rDim)  && (controlSW == 0) 
                                    k = k+1;
                                    s = rthRowNzIdx(k);
                                    [clique] = checkclique(incidenceMatrix,[p,q,r,s]);
                                    if ~isempty(clique)
                                        controlSW = 1;
                                    end
                                end
                            end
                         end
                    end
                end
            end
        end
    end
    if (p == noOfSensors-3) && (controlSW == 0)
        controlSW = -1;
    end
end

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





