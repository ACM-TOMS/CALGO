function [xMatrix,info] = SFSDPplus(sDim,noOfSensors,noOfAnchors,xMatrix0,distanceMatrix0,pars)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% A MATLAB program for executing SFSDP and showing its result
% Sunyoung Kim, Masakazu Kojima^* and Hayato Waki
% July 28, 2008
%
% Revised July 2009
% Sunyoung Kim, Masakazu Kojima^*, Hayato Waki and Makoto Yamashita
%
% * Department of Mathematical and Computing Sciences
%   Tokyo Institute of Technology
%   Oh-Okayama, Meguro, Tokyo 152-8552
%   e-mail: kojima@is.titech.ac.jp
%
% SFSDP --- A sparse version of FSDP
% Sunyoung Kim, Masakazu Kojima and Hayato Waki
% July 28, 2008
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
% distanceMatrix0 : the sparse (and noisy) distance matrix between xMatrix0(:,i) and xMatrix0(:,j);
%               distanceMatrix0(i,j) = the (noisy) distance between xMatrix0(:,i) and xMatrix0(:,j) if i < j; 
%               distanceMatrix0(i,j) = 0 if i >= j.
% pars --- parameters
%   parameters used in SeDuMi, default values: 
%       pars.eps = 1.0e-5 for sedumi;
%       pars.free = 0 for sedumi;
%       pars.fid = 0 for sedumi;
%       pars.alg, pars.theta, pars.beta, pars.vplot, pars.maxiter; 
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
%           = 2 if #anchors < sDim + 1 and no noise
%                   ---> minimizing a regularization term 
%           = 3 if #anchors < sDim + 1 and noise
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
%               constructing a FSDP relaxation when all sensors' locations are
%               given. 
%           = 0 then the original sensor network problem will be solved --- default
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Additional parameters ---> 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       pars.edgeRemoveSW
%           some different methods for eliminating edges before constructing an SDP
%           relaxation problems
%           1 --- the orignal method
%           2 --- sorting nodes according to their degrees, smaller degree nodes ---> larger degree nodes
%           3 --- sorting nodes according to their degrees, larger degree nodes ---> smaller degree nodes
%           4 --- using the MATLAB function symrcm
%           5 --- rearranging the nodes randomly
%           10 --- choosing the best one among the above methods --- default
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% <--- Additional parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% xMatrix   :   sDim x n matrix of sensors' and anchors' locations computed 
%               in the sDim dimensional space, 
%               where n is the total number of sensors and anchors, and 
%               anchors are placed in the last m_a columns 
% info      :   info from SDPsolver
%   info.cpusec     : the cpu time in seconds for the solution time.
%   info.iter       : the number of iterations. 
%   info.numerr, info.pinf, info.dinf   : see the manual of SeDuMi.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
elseif (size(xMatrix0,2) == noOfAnchors)
    if noOfAnchors < sDim+1
        fprintf('## case 1 <= noOfAnchors < sDim+1 can not be handled effectively\n'); 
        error('   Take noOfAnchors = 0 and modify distanceMatrix0.'); 
    end
    % Only anchors' location are given ---> expand the distanceMatrix0 
	% to the size sDim x (noOfSensors + noOfAnchors)
    rmsdSW = 0; 
    fprintf('## only anchor locations are given\n'); 
    xMatrix0 = [sparse(sDim,noOfSensors), xMatrix0(:,1:noOfAnchors)]; 
end 
% <--- Checking the sizes of xMatrix0
% Checking the size of distanceMatrix0 ---> 
if (size(distanceMatrix0,2) ~= noOfSensors + noOfAnchors) || (size(distanceMatrix0,1) ~= noOfSensors)
	error('## distanceMatrix0 needs to be noOfSensors x (noOfSensors + noAnchors)');
end
if norm(tril(distanceMatrix0),inf) > 1.0e-8
    fprintf('## distanceMatrix0 should be upper triangular\n');
    fprintf('   distanceMatrix0 = triu(distanceMatrix0,1)\n');
    distanceMatrix0 = triu(distanceMatrix0,1);
end
% <--- Checking the size of distanceMatrix0

if nargin < 6
    pars.analyzeData = 1; 
    pars.moreDistanceSW = 0;
else   
    if ~isfield(pars,'analyzeData')
        pars.analyzeData = 1;
    end    
    if ~isfield(pars,'moreDistanceSW')
        pars.moreDistanceSW = 0;
    end
end

% if ~isfield(pars,'SDPsolver') || isempty(pars.SDPsolver)
%      pars.SDPsolver = 'sedumi';
% end

if pars.moreDistanceSW == 1
    pars.analyzeData = 1;
end

if pars.analyzeData == 1
    [distanceMatrix0,noisyFac] = analyzeData(xMatrix0,noOfAnchors,distanceMatrix0,pars); 
    if (~isfield(pars,'noisyFac')) && (~isempty(noisyFac))
        pars.noisyFac = noisyFac;         
    end
end
if ~isfield(pars,'noisyFac')
    pars.noisyFac = 0.3;
end
noisyFac = pars.noisyFac; 

if ~isfield(pars,'eps')
    pars.eps = 1.0e-5;
end
if ~isfield(pars,'free')
    pars.free = 0;
end
if ~isfield(pars,'fid')
    pars.fid = 0;
end
% parameters for SDP relaxation ---> 
if ~isfield(pars,'minDegree')
    pars.minDegree = sDim + 2;
end
if ~isfield(pars,'objSW')
    if noisyFac < 1.0e-12
        if noOfAnchors >= sDim+1
            pars.objSW = 0;
        else
            pars.objSW = 2;
        end
    else
        if noOfAnchors >= sDim+1
            pars.objSW = 1;
        else
            pars.objSW = 3;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%
% pars.edgeRemoveSW = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%
% <--- parameters for SDP relaxation

% parameters for local refinement ---> 
if ~isfield(pars,'localMethodSW')
    pars.localMethodSW = 1;
    % = 1 --- the gradient method
    % = 0 --- no local refinement 
end
% <--- parameters for local refinement

% SDP relaxation ----> 
% 
[xMatrix,info,distanceMatrix] = SFSDP(sDim,noOfSensors,noOfAnchors,xMatrix0,distanceMatrix0,pars);
% 
% <--- SDP relaxation 

% if ~isfield(pars,'SDPsolver') || isempty(pars.SDPsolver) || strcmp(pars.SDPsolver,'sedumi')
%     SDPsolverCpuTime = info.cpusec;
%     fprintf('## cpu time for SDP solver SeDuMi = %8.2f\n',SDPsolverCpuTime);
% else
%     SDPsolverCpuTime = info.sdpaTime;
%     fprintf('## cpu time for SDP solver SDPA = %8.2f\n',SDPsolverCpuTime);
% end

fprintf('## elapsed time for SDP solver = %8.2f\n',info.SDPsolverTime);

%%%%%%%%%%%%%%%%
problemId = 10; 
%%%%%%%%%%%%%%%%

% Computing error in distance equations ---> 
[meanError,maxError] = checkDistance([xMatrix,xMatrix0(:,noOfSensors+1:noOfSensors+noOfAnchors)],distanceMatrix); 
fprintf('## mean error in dist. eq. = %6.2e, max. error in dist. eq. = %6.2e\n',meanError,maxError);
% <--- Computing error in distance equations

if rmsdSW == 1
    % Computing rmsd --->
    residualVector= reshape(xMatrix(:,1:noOfSensors) - xMatrix0(:,1:noOfSensors),1,sDim*noOfSensors);
    rmsd = norm(residualVector)/sqrt(noOfSensors);
    fprintf('## rmsd = %6.2e\n',rmsd);
    % <--- Computing rmsd
else
    rmsd = [];
end
 
% Drawing  a picture of computed and true locations of sensors ---> 
if (sDim == 2) || (sDim == 3)
    figNo=problemId*10+1; % 0+pictureNo;
    drawPicture(figNo,sDim,noisyFac,xMatrix0',[xMatrix,xMatrix0(:,noOfSensors+1:noOfSensors+noOfAnchors)]',rmsd,noOfSensors);
    fprintf('## see Figure %d\n',figNo);
else
    fprintf('\n');
end
% <--- Drawing a picture of computed and true locations of sensors

% Adujust the coordinates for anchor free cases ---> 
if ((pars.objSW == 2) || (pars.objSW == 3)) && (rmsdSW == 1) 
%    startingTime = cputime; 
    startingTime = tic; 
    xMatrix = full(xMatrix); 
	[dist,zMatrix,TRANSFORM] = procrustes(xMatrix0', xMatrix'); 
    zMatrix = zMatrix'; 
    if rmsdSW == 1
        % Computing rmsd --->
        residualVector= reshape(zMatrix(:,1:noOfSensors) - xMatrix0(:,1:noOfSensors),1,sDim*noOfSensors);
        rmsd = norm(residualVector)/sqrt(noOfSensors);
        % <--- Computing rmsd
%        fprintf('## cpu time for adjusting cordinates = %8.2f, rmsd = %6.2e',cputime - startingTime, rmsd);
        fprintf('## elapsed time for adjusting cordinates = %8.2f\n## rmsd = %6.2e',toc(startingTime), rmsd);
    else
        rmsd = [];
%        fprintf('## cpu time for adjusting cordinates = %8.2f, rmsd = %6.2e',toc(startingTime), rmsd);
        fprintf('## elapsed time for adjusting cordinates = %8.2f',cputime - startingTime);
    end
    % Drawing a picture of computed and true locations of sensors ---> 
    if (sDim == 2) || (sDim == 3)
        figNo=problemId*10+2; % 0+pictureNo;
        drawPicture(figNo,sDim,noisyFac,xMatrix0',[zMatrix,xMatrix0(:,noOfSensors+1:noOfSensors+noOfAnchors)]',rmsd,noOfSensors);
        fprintf('; see Figure %d\n',figNo);
    else
        fprintf('\n');
    end
    % <--- Drawing a picture of computed and true locations of sensors
end
% <--- Adujust the coordinates for anchor free cases 

% Refining locations of sensors by a gradient method  --->
if (pars.localMethodSW == 1) && (isempty(rmsd) || (rmsd > 1.0e-5))
    % A gradient method developed by Kim Toh 
    DD = [distanceMatrix0; distanceMatrix0(:,noOfSensors+1:noOfSensors+noOfAnchors)',sparse(noOfAnchors,noOfAnchors)];    
 %   cpuTimeStart = cputime; 
    cpuTimeStart = tic; 
    [xMatrix,Info] = ... 
        refinepositions(xMatrix(:,1:noOfSensors),xMatrix0(:,noOfSensors+1:noOfSensors+noOfAnchors),DD,3000,1.0e-8); 
%    fprintf('## cpu time for a gradient method = %8.2f\n',cputime - cpuTimeStart);
    fprintf('## elapsed time for a gradient method = %8.2f\n',toc(cpuTimeStart));
    xMatrix = [xMatrix,xMatrix0(:,noOfSensors+1:noOfSensors+noOfAnchors)];
    % Computing error in distance equations ---> 
    [meanError,maxError] = checkDistance(xMatrix,distanceMatrix); 
    % <--- Computing error in distance equations
    fprintf('## mean error in dist. eq. = %6.2e, max. error in dist. eq. = %6.2e\n',meanError,maxError);    
    % Adujust the coordinates for anchor free cases ---> 
    if (pars.objSW == 2) || (pars.objSW == 3) && (rmsdSW == 1) 
        % cpuTimeStart = cputime; 
        cpuTimeStart = tic; 
        xMatrix = full(xMatrix);
        [dist,zMatrix] = procrustes(xMatrix0', xMatrix');
        xMatrix = zMatrix';
%        fprintf('## cpu time for adjusting cordinates = %8.2f\n',cputime-cpuTimeStart);
        fprintf('## elapsed time for adjusting cordinates = %8.2f\n',toc(cpuTimeStart));
    % <--- Adujust the coordinates for anchor free cases 
        % Computing rmsd --->
    end
    if (rmsdSW == 1)
        residualVector= reshape(xMatrix(:,1:noOfSensors) - xMatrix0(:,1:noOfSensors),1,sDim*noOfSensors);
        rmsd2 = norm(residualVector)/sqrt(noOfSensors);
        % <--- Computing rmsd
        fprintf('## rmsd = %6.2e\n',rmsd2);
    end
    % <--- Computing rmsd
    % Drawing a picture of computed and true locations of sensors ---> 
    if (sDim == 2) || (sDim == 3)
        figNo=problemId*10+3; % 0+pictureNo;
        drawPicture(figNo,sDim,noisyFac,xMatrix0',xMatrix',rmsd,noOfSensors);
        fprintf('## see Figure %d\n',figNo);
    end
    % <--- Drawing a picture of computed and true locations of sensors    
else
    rmsd2 = [];
end
% <--- Refining locations of sensors by a gradient method

% fprintf('\n'); 

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [distanceMatrix,standardDeviation] = analyzeData(xMatrix0,noOfAnchors,distanceMatrix,pars)

standardDeviation = [];

if nargin < 4
    pars.moreDistanceSW = 0;
elseif ~isfield(pars,'moreDistanceSW')
    pars.moreDistanceSW = 0;
end
    
sDim = size(xMatrix0,1);
noOfSensors = size(xMatrix0,2) - noOfAnchors;
distanceMatrix = triu(distanceMatrix,1);
distanceMatrix = distanceMatrix(1:noOfSensors,:); 
noOfEdges = nnz(distanceMatrix); 
noOfEdges2 = nnz(distanceMatrix(:,1:noOfSensors)); 

degreeVector = sum(spones([distanceMatrix(:,1:noOfSensors)+distanceMatrix(:,1:noOfSensors)',...
    distanceMatrix(:,noOfSensors+1:noOfSensors+noOfAnchors)]),2); 
minDeg = full(min(degreeVector')); 
maxDeg = full(max(degreeVector')); 
averageDeg = full(sum(degreeVector')/noOfSensors); 

fprintf('## sDim = %d, noOfSensors = %d, noOfAnchors = %d\n',sDim,noOfSensors,noOfAnchors);
fprintf('## the number of dist. eq. between two sensors  = %d\n',noOfEdges2);
fprintf('## the number of dist. eq. between a sensor & an anchor = %d\n',noOfEdges-noOfEdges2);
fprintf('## the min., max. and ave. degrees over sensor nodes = %d, %d, %6.2f\n',minDeg,maxDeg,averageDeg);

if norm(xMatrix0(:,1:noOfSensors),inf) > 0
    sensorDataSW = 1;
else
    sensorDataSW = 0;
end
    
if sensorDataSW == 1
    coordinateRange = zeros(sDim,2);
    for i=1:sDim
        coordinateRange(i,1) = min(xMatrix0(i,:));
        coordinateRange(i,2) = max(xMatrix0(i,:));
    end
    for i=1:sDim
        fprintf('## %+9.4e <= x(%d) <= %+9.4e\n',coordinateRange(i,1),i,coordinateRange(i,2));
    end
else
    fprintf('## no location for sensors is given\n'); 
end

% minRadioRange = 1.0e10;
maxRadioRange = 0.0;
% maxErrorRatio = 0.0;
% sumErrorRatio = 0.0;
noOfEdges = 0; 
totalSquaredErrorRatio = 0.0; 

if sensorDataSW == 1
    for p=1:noOfSensors
        nzIdx = find(distanceMatrix(p,:) > 0);
        noOfEdges = noOfEdges + length(nzIdx);
        for q=nzIdx
            trueDist = norm(xMatrix0(:,p)-xMatrix0(:,q));
            measuredDist = distanceMatrix(p,q);
            %        minRadioRange = min(minRadioRange,trueDist);
            maxRadioRange = max(maxRadioRange,trueDist);
            errorRatio = abs(measuredDist-trueDist)/trueDist;
            %        maxErrorRatio = max(maxErrorRatio,errorRatio);
            %        sumErrorRatio = sumErrorRatio + errorRatio;
            totalSquaredErrorRatio = totalSquaredErrorRatio + errorRatio^2;
        end
    end
    standardDeviation = sqrt(totalSquaredErrorRatio/(noOfEdges-1));
    fprintf('## the max. radio range = %9.4e, the estimated noisy factor = %9.4e\n',maxRadioRange,standardDeviation);

    % minRadioRange
    % maxRadioRange
    % maxErrorRatio
    % noOfEdges
    % averageErrorRatio = sumErrorRatio/noOfEdges
    % standardDeviation = sqrt(totalSquaredErrorRatio/(noOfEdges-1))
    if pars.moreDistanceSW == 1
        fprintf('## all other edges (with noise) are added within the max. radio range = %9.4e:\n',maxRadioRange);
        randSeed = 3201;
        randn('seed',randSeed);
        radioRange = maxRadioRange;
        noisyFac = standardDeviation;
        totalSquaredErrorRatio = 0.0;
        for p=1:noOfSensors
            for q=p+1:noOfSensors+noOfAnchors
                trueDist = norm(xMatrix0(:,p)-xMatrix0(:,q));
                if distanceMatrix(p,q) == 0
                    if trueDist <= radioRange
                        epsilon0 = randn(1,1)*noisyFac;
                        rate = max([1+epsilon0,0.1]);
                        distanceMatrix(p,q) = rate*trueDist;
                        %            distanceMatrix(p,q) = trueDist;
                        totalSquaredErrorRatio = totalSquaredErrorRatio + epsilon0^2;
                    end
                else
                    measuredDist = distanceMatrix(p,q);
                    errorRatio = abs(measuredDist-trueDist)/trueDist;
                    totalSquaredErrorRatio = totalSquaredErrorRatio + errorRatio^2;
                end
            end
        end
        noOfEdges = nnz(distanceMatrix);
        standardDeviation = sqrt(totalSquaredErrorRatio/(noOfEdges-1));
        noOfEdges = nnz(distanceMatrix);
        noOfEdges2 = nnz(distanceMatrix(:,1:noOfSensors));

        degreeVector = sum(spones([distanceMatrix(:,1:noOfSensors)+distanceMatrix(:,1:noOfSensors)',...
            distanceMatrix(:,noOfSensors+1:noOfSensors+noOfAnchors)]),2);
        minDeg = full(min(degreeVector'));
        maxDeg = full(max(degreeVector'));
        averageDeg = full(sum(degreeVector')/noOfSensors);

        fprintf('   the computed noisy factor= %9.4e\n',standardDeviation);
        fprintf('   the number of dist. eq. between two sensors  = %d\n',noOfEdges2);
        fprintf('   the number of dist. eq. between a sensor & an anchor = %d\n',noOfEdges-noOfEdges2);
        fprintf('   the min., max. and ave. degrees over sensor nodes = %d, %d, %6.2f\n',minDeg,maxDeg,averageDeg);
    end
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [meanError,maxError] = checkDistance(xMatrix,distanceMatrix)
noOfSEnsors = size(xMatrix,1);
totalError = 0.0; 
maxError = 0.0;
for p=1:noOfSEnsors
    nzIdx = find(distanceMatrix(p,:)); 
    for q = nzIdx
        error = abs(norm(xMatrix(:,p)-xMatrix(:,q)) - distanceMatrix(p,q)); 
        totalError = totalError + error;
        maxError = max(maxError,error);
    end
end
meanError = totalError/nnz(distanceMatrix);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function drawPicture(figNo,sDim,noisyFac,xData,xComputed,rmsd,noOfSensors)
if sDim == 2
    draw2dPicture(figNo,noisyFac,xData,xComputed,rmsd,noOfSensors)
elseif sDim == 3
    draw3dPicture(figNo,noisyFac,xData,xComputed,rmsd,noOfSensors); 
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function draw2dPicture(figNo,noisyFac,xData,xComputed,rmsd,noOfSensors)

epsilon = 11.0e-7;
tol = max([noisyFac,epsilon]); 

if isempty(rmsd)
    xMin = min(xComputed(:,1)');
    xMax = max(xComputed(:,1)');
    yMin = min(xComputed(:,2)');
    yMax = max(xComputed(:,2)');
    len = size(xData,1);
    figure(figNo);
    %pause
    %axis([-0.6 0.6 -0.6 0.6]);
    %axis([0.0 1.0 0.0 1.0]);
    
%    [xMin xMax yMin yMax]
    
    axis([xMin xMax yMin yMax]);

    plot(xComputed(1:noOfSensors,1),xComputed(1:noOfSensors,2),'r*');
    hold on;
    plot(xData(noOfSensors+1:len,1),xData(noOfSensors+1:len,2),'bD');
    if (mod(figNo,10) == 1)
        text(xMin+0.2, yMin-0.08*(yMax-yMin),'* : the sensor locations computed by SFSDP');
    elseif (mod(figNo,10) == 2)
        text(xMin-0.05, yMin-0.08*(yMax-yMin),'* : the sensor locations computed by SFSDP + the procrustes funct. from K. Toh');
    elseif (mod(figNo,10) == 3)
        text(xMin-0.05, yMin-0.08*(yMax-yMin),'* : the sensor locations computed by SFSDP + the refinepositions funct. by K. Toh');
    end
    hold off
else
    xMin = min(xData(:,1)');
    xMax = max(xData(:,1)');
    yMin = min(xData(:,2)');
    yMax = max(xData(:,2)');
    % select points of errors
    xOri = [];
    xFou = [];
    for i=1:noOfSensors
        if norm(xData(i,:)-xComputed(i,:)) > 0.002
            xOri =[xOri;xData(i,:)];
            xFou =[xFou;xComputed(i,:)];
        end
    end

    len = size(xData,1);

    figure(figNo);
    
    plot(xData(1:noOfSensors,1),xData(1:noOfSensors,2),'gO');

    hold on;
    %pause
    %axis([-0.6 0.6 -0.6 0.6]);
    %axis([0.0 1.0 0.0 1.0]);
    axis([xMin xMax yMin yMax]);

    plot(xComputed(1:noOfSensors,1),xComputed(1:noOfSensors,2),'r*');

    %
    plot(xData(noOfSensors+1:len,1),xData(noOfSensors+1:len,2),'bD');

    %title('text',' sensors ',noOfSensors, ' anchors ',len-noOfSensors, 'noise = ',tol);
    % rmsd

    
    for i=1:size(xOri,1)
        xx = [xOri(i,1),xFou(i,1)];
        yy = [xOri(i,2),xFou(i,2)];
        plot(xx,yy,'b-');
    end

    % annotation('textbox',[0.2 0.02 0.6 0.05]);
    if (mod(figNo,10) == 1)
        text(xMin+0.2, yMin-0.08*(yMax-yMin),'O : Sensor true locations   vs   * : the ones computed by SFSDP');
    elseif (mod(figNo,10) == 2)
        text(xMin-0.05, yMin-0.08*(yMax-yMin),'O : Sensor true locations   vs   * : the ones computed by SFSDP + the procrustes funct. from K. Toh');
    elseif (mod(figNo,10) == 3)
        text(xMin-0.05, yMin-0.08*(yMax-yMin),'O : Sensor true locations   vs   * : the ones computed by SFSDP + the refinepositions funct. by K. Toh');
    end
    hold off
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function draw3dPicture(figNo,noisyFac,xData,xComputed,rmsd,noOfSensors); 

epsilon = 11.0e-7;
tol = max([noisyFac,epsilon]); 

if ~isempty(rmsd)
    % select points of errors
    xOri = [];
    xFou = [];
    for i=1:noOfSensors
        if norm(xData(i,:)-xComputed(i,:)) > 0.002
            xOri =[xOri;xData(i,:)];
            xFou =[xFou;xComputed(i,:)];
        end
    end

    len = size(xData,1);
    figure(figNo);
    plot3(xData(1:noOfSensors,1),xData(1:noOfSensors,2),xData(1:noOfSensors,3),'gO');
    hold on;
    %pause
    % x1 = -0.6:0.02:0.6;
    % y1 = -0.6:0.02:0.6;
    % zn = length(x1);
    % z1 = -0.6*ones(zn,zn);
    %z1(zn,zn) = 0.6;
    % [X,Y] = meshgrid(x1,y1);
    % mesh(X,Y,z1);

    % xMin = round(min(xData(:,1)'));
    % xMax = round(max(xData(:,1)'));
    % yMin = round(min(xData(:,2)'));
    % yMax = round(max(xData(:,2)'));
    xMin = min(xData(:,1)');
    xMax = max(xData(:,1)');
    yMin = min(xData(:,2)');
    yMax = max(xData(:,2)');
    zMin = min(xData(:,3)');
    zMax = max(xData(:,3)');

    %axis([-0.6 0.6 -0.6 0.6 -0.6 0.6]);
    % axis([0.0 1.0 0.0 1.0 0.0 1.0]);
    axis([xMin xMax yMin yMax zMin zMax]);

    plot3(xComputed(1:noOfSensors,1),xComputed(1:noOfSensors,2),xComputed(1:noOfSensors,3),'r*');

    plot3(xData(noOfSensors+1:len,1),xData(noOfSensors+1:len,2),xData(noOfSensors+1:len,3),'bD');

    grid on
    %title('text',' sensors ',noOfSensors, ' anchors ',len-noOfSensors, 'noise = ',tol);
    % rmsd
    for i=1:size(xOri,1)
        xx = [xOri(i,1),xFou(i,1)];
        yy = [xOri(i,2),xFou(i,2)];
        zz = [xOri(i,3),xFou(i,3)];
        plot3(xx,yy,zz,'b-');
    end
    % annotation('textbox',[0.2 0.02 0.6 0.05]);
    % if (mod(figNo,10) == 1)
    %     text(-0.6, 0.10,'O : Sensor true locations   vs   * : the ones computed by SFSDP');
    % elseif (mod(figNo,10) == 2)
    %     text(-0.8, 0.22,'O : Sensor true locations   vs   * : the ones computed by SFSDP + the procrustes funct. from K. Toh');
    % elseif (mod(figNo,10) == 3)
    %     text(-0.8, 0.22,'O : Sensor true locations   vs   * : the ones computed by SFSDP + the refinepositions funct. by K. Toh');
    % end
    hold off
else
    len = size(xData,1);
    figure(figNo);
    xMin = min(xComputed(:,1)');
    xMax = max(xComputed(:,1)');
    yMin = min(xComputed(:,2)');
    yMax = max(xComputed(:,2)');
    zMin = min(xComputed(:,3)');
    zMax = max(xComputed(:,3)');

    %axis([-0.6 0.6 -0.6 0.6 -0.6 0.6]);
    % axis([0.0 1.0 0.0 1.0 0.0 1.0]);
    axis([xMin xMax yMin yMax zMin zMax]);

    plot3(xComputed(1:noOfSensors,1),xComputed(1:noOfSensors,2),xComputed(1:noOfSensors,3),'r*');
    hold on;
    plot3(xData(noOfSensors+1:len,1),xData(noOfSensors+1:len,2),xData(noOfSensors+1:len,3),'bD');
    grid on
    
    % if (mod(figNo,10) == 1)
    %     text(-0.6, 0.10,'O : Sensor true locations   vs   * : the ones computed by SFSDP');
    % elseif (mod(figNo,10) == 2)
    %     text(-0.8, 0.22,'O : Sensor true locations   vs   * : the ones computed by SFSDP + the procrustes funct. from K. Toh');
    % elseif (mod(figNo,10) == 3)
    %     text(-0.8, 0.22,'O : Sensor true locations   vs   * : the ones computed by SFSDP + the refinepositions funct. by K. Toh');
    % end
    hold off
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [d, Z, transform] = procrustes(X, Y)
%
% Thanks to Kim Chuan Toh
% 
%PROCRUSTES Procrustes Analysis
%   D = PROCRUSTES(X, Y) determines a linear transformation (translation,
%   reflection, orthogonal rotation, and scaling) of the points in the
%   matrix Y to best conform them to the points in the matrix X.  The
%   "goodness-of-fit" criterion is the sum of squared errors.  PROCRUSTES
%   returns the minimized value of this dissimilarity measure in D.  D is
%   standardized by a measure of the scale of X, given by
%
%      sum(sum((X - repmat(mean(X,1), size(X,1), 1)).^2, 1))
%
%   i.e., the sum of squared elements of a centered version of X.  However,
%   if X comprises repetitions of the same point, the sum of squared errors
%   is not standardized.
%
%   X and Y are assumed to have the same number of points (rows), and
%   PROCRUSTES matches the i'th point in Y to the i'th point in X.  Points
%   in Y can have smaller dimension (number of columns) than those in X.
%   In this case, PROCRUSTES adds columns of zeros to Y as necessary.
%
%   [D, Z] = PROCRUSTES(X, Y) also returns the transformed Y values.
%
%   [D, Z, TRANSFORM] = PROCRUSTES(X, Y) also returns the transformation
%   that maps Y to Z.  TRANSFORM is a structure with fields:
%      c:  the translation component
%      T:  the orthogonal rotation and reflection component
%      b:  the scale component
%   That is, Z = TRANSFORM.b * Y * TRANSFORM.T + TRANSFORM.c.
%
%   Examples:
%
%      % Create some random points in two dimensions
%      X = normrnd(0, 1, [10 2]);
%
%      % Those same points, rotated, scaled, translated, plus some noise
%      S = [0.5 -sqrt(3)/2; sqrt(3)/2 0.5]; % rotate 60 degrees
%      Y = normrnd(0.5*X*S + 2, 0.05, size(X));
%
%      % Conform Y to X, plot original X and Y, and transformed Y
%      [d, Z, tr] = procrustes(X,Y);
%      plot(X(:,1),X(:,2),'rx', Y(:,1),Y(:,2),'b.', Z(:,1),Z(:,2),'bx');
%
%   See also FACTORAN, CMDSCALE.

%   References:
%     [1] Seber, G.A.F., Multivariate Observations, Wiley, New York, 1984.
%     [2] Bulfinch, T., The Age of Fable; or, Stories of Gods and Heroes,
%         Sanborn, Carter, and Bazin, Boston, 1855.

%   Copyright 1993-2002 The MathWorks, Inc.
%   $Revision: 1.2 $  $Date: 2002/01/17 21:31:45 $

[n, m]   = size(X);
[ny, my] = size(Y);

if ny ~= n
    error('X and Y must have the same number of rows (points).');
elseif my > m
    error('Y cannot have more columns (variables) than X.');
end

% center at the origin
muX = mean(X,1);
muY = mean(Y,1);
X0 = X - repmat(muX, n, 1);
Y0 = Y - repmat(muY, n, 1);

ssqX = sum(X0.^2,1);
ssqY = sum(Y0.^2,1);
constX = all(ssqX <= abs(eps*n*muX).^2);
constY = all(ssqY <= abs(eps*n*muY).^2);

if ~constX & ~constY
    % the "centered" Frobenius norm
    normX = sqrt(sum(ssqX)); % == sqrt(trace(X0*X0'))
    normY = sqrt(sum(ssqY)); % == sqrt(trace(Y0*Y0'))

    % scale to equal (unit) norm
    X0 = X0 / normX;
    Y0 = Y0 / normY;

    % make sure they're in the same dimension space
    if my < m
        Y0 = [Y0 zeros(n, m-my)];
    end

    % optimum rotation matrix of Y
    A = X0' * Y0;
    [L, D, M] = svd(A);
    T = M * L';

    % optimum (symmetric in X and Y) scaling of Y
    trsqrtAA = sum(diag(D)); % == trace(sqrtm(A'*A))

    % the standardized distance between X and bYT+c
    d = 1 - trsqrtAA.^2;
    if nargout > 1
        Z = normX*trsqrtAA * Y0 * T + repmat(muX, n, 1);
    end
    if nargout > 2
        if my < m
            T = T(1:my,:);
        end
        b = trsqrtAA * normX / normY;
        transform = struct('T',T, 'b',b, 'c',repmat(muX - b*muY*T, n, 1));
    end

% the degenerate cases: X all the same, and Y all the same
elseif constX
    d = 0;
    Z = repmat(muX, n, 1);
    T = eye(my,m);
    transform = struct('T',T, 'b',0, 'c',Z);
else % ~constX & constY
    d = 1;
    Z = repmat(muX, n, 1);
    T = eye(my,m);
    transform = struct('T',T, 'b',0, 'c',Z);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Thanks to Kim Chuan Toh
% 
%%*************************************************************************
%% refinepositions: steepest descent with back-tracking line search to 
%%                  minimize 
%%
%% f(X) = sum_{j<=k: djk is given}  (norm(xj-xk)-djk)^2 
%%        + sum_{j,k: djk is given} (norm(xj-ak)-djk)^2
%%
%% input: X0     = [sensor position in column format]
%%        anchor = [anchor position in column format] 
%%        DD     = [sensor-sensor distance, senor-anchor distance
%%                  anchor-sensor distance,  0]
%% (optional) maxit = maximum number of gradient iterations allowed 
%% (optional) tol   = stopping tolerance  
%%
%% output: Xiter = [refined sensor position]
%%         Info.objective  = history of objective values
%%         Info.gradnorm   = history of norm of the gradients
%%         Info.steplength = history of step-lengths
%%         Info.cputime    = cputime taken
%% 
%% child functions: gradfun2, objfun2
%%*************************************************************************

  function [Xiter,Info] = refinepositions(X0,anchor,DD,maxit,tol)

  if (nargin < 4); maxit = 3000; end
  if (nargin < 5); tol = 1e-9; end

  ttime = cputime; 
 
  [dummy,n] = size(X0);
  [dummy,m] = size(anchor); 
  if (size(DD,2) ~= n+m) 
     error('gradescent: dimension of X0 or DD not correct')
  end
  D1 = DD(1:n,1:n); 
  [II,JJ,dd] = find(triu(D1)); 
  if (size(dd,1) < size(dd,2)); dd = dd'; end
  ne = length(II); 
  S  = sparse(II,1:ne,1,n+m,ne)-sparse(JJ,1:ne,1,n+m,ne);   
  if (m > 0)
    %
    % 2008-06-13 Waki
    % change [1:m] into (1:m)
    %
     D2 = DD(1:n,n+(1:m)); 
     [I2,J2,d2] = find(D2); 
     if (size(d2,1) < size(d2,2)); d2 = d2'; end
     ne2 = length(I2);
     S2 = sparse(I2,1:ne2,1,n+m,ne2)-sparse(n+J2,1:ne2,1,n+m,ne2);
     S  = [S,S2]; 
     dd = [dd; d2]; 
  end
  Sg = S'; Sg = [Sg(1:length(dd),1:n), sparse(length(dd),m)]; 
  X0 = [X0, anchor]; 
  obj = objfun2(X0,S,dd); 
  gradx = gradfun2(X0,S,dd,Sg);
%%
  Info.objective  = zeros(1,maxit); 
  Info.gradnorm   = zeros(1,maxit); 
  Info.steplength = zeros(1,maxit); 
  Info.objective(1)  = obj; 
  Info.gradnorm(1)   = sqrt(max(sum(gradx.*gradx)));
%%
  Xiter = X0; objold = obj; 
  for iter = 1:maxit
     objnew = inf;
     gradx = gradfun2(Xiter,S,dd,Sg);
     alpha = 0.2; count = 0;
    %
    % 2008-06-13 Waki
    % change & into &&
    %
     while (objnew > objold) && (count < 20)
        alpha = 0.5*alpha; count = count + 1; 
        Xnew  = Xiter - alpha*gradx;
        objnew = objfun2(Xnew,S,dd);
     end 
     Xiter = Xnew;
     Info.objective(iter+1)  = objnew;
     Info.gradnorm(iter+1)   = sqrt(max(sum(gradx.*gradx)));
     Info.steplength(iter+1) = alpha; 
     if (abs(objnew-objold)/(1+abs(objold)) < tol); break; end 
     objold = objnew; 
  end
  if (m > 0); Xiter = Xiter(:,1:n); end
  ttime = cputime-ttime; 
  Info.cputime = ttime; 

%%*************************************************************************
%% Find the function value
%% f(X) = sum_{j<=k} (norm(xj-xk)-djk)^2 + sum_{j,k} (norm(xj-ak)-djk)^2
%%
%% input: X = [sensor position, anchor]
%%*************************************************************************
    %
    % 2008-06-13 Waki
    % Remove colon
    %
  function  objval = objfun2(X,S,dd)

  Xij = X*S; 
  normXij = sqrt(sum(Xij.*Xij))';  
  objval = norm(normXij-dd)^2;

%%*************************************************************************
%% Find the gradient of the function 
%% f(X) = sum_{j<=k} (norm(xj-xk)-djk)^2 + sum_{j,k} (norm(xj-ak)-djk)^2
%%
%% input: X = [sensor position, anchor]
%%*************************************************************************
    %
    % 2008-06-13 Waki
    % Remove colon
    %
  function  G = gradfun2(X,S,dd,Sg)

  ne = length(dd);
  Xij = X*S;
  normXij = sqrt(sum(Xij.*Xij))'+eps;
  tmp = 1-dd./normXij;
  G = Xij*spdiags(2*tmp,0,ne,ne);
  G = G*Sg;
%%*************************************************************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    
    

