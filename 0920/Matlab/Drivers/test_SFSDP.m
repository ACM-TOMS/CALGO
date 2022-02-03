function test_SFSDP(sDim,noisyFac,radiorange,noOfSensors,anchorType,noOfAnchors,randSeed);
%%
% To use this program, SFSDP.m, SFDPplus.m, generateProblem.m and sdpa.7.3.1 (or sedumi) 
% are necessary. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% A MATLAB program for testing SFSDP
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
% P. Biswas and Y. Ye (2004) gSemidefinite programming for ad hoc wireless 
% sensor network localization,h in Proceedings of the third international 
% symposium on information processing in sensor networks, ACM press, 46-54.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Acknowledgments:
% 
% The authors are grateful to Professor Yinyu Ye for the original 
% version of FSDP, and Professor Kim Chuan Toh 
% for MATLAB programs refineposition.m and procrustes.m, 
% and valuable suggestions. 
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
% Input 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sDim :    the dimension of the space in which sensors and anchors are
%           located; sDim is either 2 or 3. 
%           When sDim = 2, sensors and anchors are located in [0,1]^2.
%           When sDim = 3, sensors and anchors are located in [0,1]^3.
% noisyFac :    noisy factor 
%               = 0 --- no noise 
%               = \sigma > 0 --- noise with the distribution N(0,\sigma).
% radiorange :  radio range; 
%               If \|x_p - \x_q\| <= radio range, a distance (with noise) is given between x_p and x_q. 
%               If \|x_p - \x_q\| > radio range, no distance is given between x_p and x_q. 
%               Here all sensors are placed in [0,1]^{sDim} randomly.                
% noOfSensors : the number of sensors.
% anchorType :	parameters on how we locate anchors. 
%   anchorType = 0
%       all anchors are placed at grid points on the boundary and interior
%       of [0,1]^{sDim}
%   anchorType = 1
%       all anchors are placed at grid points in the interior of
%       [0,1]^{sDim}
% 	anchorType = 2
%       all anchors are placed randomly in [0,1]^{sDim}
%   anchorType = 3
%       sDim + 1 anchors on the origin and the coordinate axis
%   anchorType = 4
%       sDim + 1 anchors near the center
%   anchorType = 10
%       no anchor
% noOfAnchors : the number of anchors
% randSeed  :   a random seed number for a random distribution of sensors and
%               anchors when anchorType = 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% rmsd : the root mean square distance of the solution computed by SFSDP. 
% rmsd2 : the root mean square distance of the solution computed by SFSDP + the gradient method.
% sedumiCpuTime : the SeDuMi cpu time. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
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
%             bound the error epsilon_{ij}^+ and epsilon_{ij}^-. 
%       pars.localMethodSW 
%           = 1 --->    apply the gradient method to improve the solution
%                       obtained by SFSDP. 
%           = 0 --->    No application of the gradient method
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Choose an SDP solver,  sedumi or sdpa 
%  
% pars.SDPsolver = 'sedumi';
pars.SDPsolver = 'sdpa'; 
%

pars.sparseSW = 1;

% Checking anchorType and noOfAnchors  --->
if ((anchorType < 0) || (anchorType > 4)) && (anchorType ~= 10)
    error('anchorType should be 0, 1, 2, 3, 4 and 10.');
elseif mod(anchorType,1) > 0
    error('anchorType should be 0, 1, 2, 3, 4 and 10.');
elseif (anchorType == 0) || (anchorType == 1)
    if noOfAnchors < power(2,sDim)
        error('When anchorType = 0 or 1, noOfAnchors has to be n^{sDim} for some integer n >= 2.'); 
    else
        gridForAnchors = power(noOfAnchors,1/sDim); 
        if mod(gridForAnchors,1) > 0
            error('When anchorType = 0 or 1, noOfAnchors has to be n^{sDim} for some integer n >= 2.'); 
        end
    end    
elseif (anchorType == 2)
    if noOfAnchors < sDim+1
       error('noOfAnchors should be not less than sDim+1 when anchorType = 2.');
    end
elseif (anchorType == 3) || (anchorType == 4)
    if noOfAnchors ~= sDim+1
       error('noOfAnchors should be sDim+1 when anchorType = 3 or 4.');
    end 
elseif anchorType == 10
    if noOfAnchors > 0
        error('noOfAnchors should be 0 when anchorType = 10.');
    end
end
% <--- Checking anchorType and noOfAnchors

% if nargin <= 7
%     problemId = 0;
% end

problemId = 10;
 
% All sensors and anchors are located in the [0,1] x [0,1] or
% [0,1] x [0,1] x [0,1] region. 

% if problemId >= 1 
%     fprintf('## problemId = %d\n',problemId); 
% end

% Generating a sensor network problem ---> 
% 

% startingTime = cputime; 

[xMatrix0,distanceMatrix] ...
	= generateProblem(sDim,noisyFac,radiorange,noOfSensors,anchorType,noOfAnchors,randSeed);

pars.noisyFac = noisyFac; 

noOfEdges = nnz(distanceMatrix); 
noOfEdges2 = nnz(distanceMatrix(:,1:noOfSensors)); 
degreeVector = sum(spones([distanceMatrix(:,1:noOfSensors)+distanceMatrix(:,1:noOfSensors)',...
    distanceMatrix(:,noOfSensors+1:noOfSensors+noOfAnchors)]),2); 
minDeg = full(min(degreeVector')); 
maxDeg = full(max(degreeVector')); 
averageDeg = full(sum(degreeVector')/noOfSensors); 

fprintf('## sDim = %d, noOfSensors = %d, anchorType = %d, noOfAnchors = %d\n',sDim,noOfSensors,anchorType,noOfAnchors);
fprintf('## radiorange = %6.2e, noisyFac = %6.2e, randSeed = %d\n',radiorange,noisyFac,randSeed);
fprintf('## the number of dist. eq. between two sensors  = %d\n',noOfEdges2);
fprintf('## the number of dist. eq. between a sensor & an anchor = %d\n',noOfEdges-noOfEdges2);
fprintf('## the min., max. and ave. degrees over sensor nodes = %d, %d, %6.2f\n',minDeg,maxDeg,averageDeg);

% parameters for SDPsolver, sedumi or sdpa ---> 
%     if strcmp(pars.SDPsolver,'sedumi')
%         pars.eps = 1.0e-5; 
%     else
    pars.eps = 1.0e-5;
%     end
    pars.free = 0;
    pars.fid = 0;
% <--- parameters for SeDuMi 
% parameters for SDP relaxation ---> 
    pars.minDegree = sDim + 2; % Increase pars.minDegree for a better accuracy
    if noisyFac < 1.0e-12
        if noOfAnchors >= sDim+1
            pars.objSW = 0;
        else
            pars.objSW = 2;
        end   
    %    pars.objSW = 2;
    else
        if noOfAnchors >= sDim+1
            pars.objSW = 1;
        else
            pars.objSW = 3;
    %        pars.minDegree = sDim + 2; 
        end   
    %    pars.objSW = 3;
    end
    if noisyFac > 0
        pars.noisyFac = noisyFac; 
    end    
% <--- parameters for SDP relaxation
% parameters for local refinement ---> 
    pars.localMethodSW = 1;
    % = 1 --- the gradient method
    % = 0 --- no local refinement 
% <--- parameters for local refinement

pars.analyzeData = 0; 
pars.moreDistanceSW = 0;

%%%%%%
% pars.minDegree = sDim+2; 
%%%%%%

[xMatrix,info] = SFSDPplus(sDim,noOfSensors,noOfAnchors,xMatrix0,distanceMatrix,pars); 

return

