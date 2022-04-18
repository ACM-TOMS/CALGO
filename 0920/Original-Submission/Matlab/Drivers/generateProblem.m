function [xMatrix0,distanceMatrix0] = generateProblem(sDim,noisyFac,... 
    radiorange,noOfSensors,anchorType,noOfAnchors,randSeed)
% A MATLAB program for generating a sensor network localization problem
% Sunyoung Kim, Masakazu Kojima^* and Hayato Waki
% July 28, 2008
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%               Here all sensors are placed in [0,1]^{sDim} randomely.                
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% xMatrix0 :    an sDim x (noOfSensors + noOfAnchors) matrix to represent 
%               locations of sensors and anchors, where noOfSensors sensors
%               are placed in the first noOfSensors columns, and
%               noOfAnchors are place in the last noOfAnchors columns. 
% distanceMatrix0 :  a noOfSensors x (noOfSensors + noOfAnchors) upper
%                   triangular matrix to represent distances from sensors
%                   to sensors and anchors, where the (p,q)th element
%                   d_{pq} denotes the distance from the pth sensor to the qth
%                   sensor or anchor. When noisyFac = 0,
%                       d_{pq} = norm(xMatrix0(:,p) - xMatrix0(:,q))
%                           if norm(xMatrix0(:,p) - xMatrix0(:,q)) <= radiorange,
%                       d_{pq} = 0 
%                           if norm(xMatrix0(:,p) - xMatrix0(:,q)) > radiorange,
%                   When \sigma = noisyFac > 0,
%                       d_{pq} = max(1+epsilon,0.1)*(norm(xMatrix0(:,p) - xMatrix0(:,q))
%                           if norm(xMatrix0(:,p) - xMatrix0(:,q)) <= radiorange,
%                       d_{pq} = 0 
%                           if norm(xMatrix0(:,p) - xMatrix0(:,q)) > radiorange,
%                   Here epsilon is a random number chosen from N(0,\sigma). 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Other parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% degreeLB1 :   if a sensor p is already connected to at least degreeLB1 anchors 
%               then no edge from the sensor p to any other anchors is
%               added. 
% degreeLB1 = 1.0e10; i.e., degreeLB1 = \infty, all edges between a sensor and 
%               an anchors whose distance is not greater than the radiorange are 
%               generated. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% degreeLB2 :    if degree(p) \geq degreeLB and degree(q) \beq degreeLB hold at
%                sensors p and q, then the edge is not added between them. 
% degreeLB2 = 1.0e10;  % i.e., degreeLB = \infty, all edges whose distance are 
%                     % not greater than the radiorange are generated. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if noisyFac == 0
    degreeLB1 = sDim+1; 
    degreeLB2 = (sDim+1) * 5; 
else 
    degreeLB1 = (sDim+1)*2; 
    degreeLB2 = (sDim+1) * 6; 
end

xbd = 1.0;

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

startingTime = tic; % cputime; 

rand('seed',randSeed);

% placing anchors on grids ---> 
if (anchorType == 0) || (anchorType == 1)
    gridForAnchors = power(noOfAnchors,1/sDim); 
    if noOfAnchors == 1
        noOfAnchors = 1;
        anchorMatrix = ones(sDim,1)*0.5;
        xMatrix0 = [xbd*rand(sDim,noOfSensors),anchorMatrix];
    else % noOfAnchors > 1
        xMatrix0 = xbd*rand(sDim,noOfSensors+noOfAnchors);
        anchorMatrix = sparse(sDim,noOfAnchors);
        h0 = xbd/(gridForAnchors+2);
        if anchorType == 0
            h = xbd / (gridForAnchors - 1);
            origin = zeros(sDim,1);
        else
            h = (xbd-h0) / (gridForAnchors-1);
            origin = (h0/2) * ones(sDim,1);
        end
        pointer = 0;
        if sDim == 2
            for i = 1:gridForAnchors
                for j=1:gridForAnchors
                    pointer = pointer+1;
                    anchorMatrix(1,pointer) = (i-1)*h + origin(1,1);
                    anchorMatrix(2,pointer) = (j-1)*h + origin(2,1);
                    % fprintf('%s Anchor %d = (%17.15f,%17.15f)\n',char(37),pointer,anchorMatrix(1,pointer),anchorMatrix(2,pointer));
                end
            end
        elseif sDim == 3
            for i = 1:gridForAnchors
                for j=1:gridForAnchors
                    for k=1:gridForAnchors
                        pointer = pointer+1;
                        anchorMatrix(1,pointer) = (i-1)*h + origin(1,1);
                        anchorMatrix(2,pointer) = (j-1)*h + origin(2,1);
                        anchorMatrix(3,pointer) = (k-1)*h + origin(3,1);
                        % fprintf('%s Anchor %d = (%17.15f,%17.15f,%17.15f)\n',char(37),pointer,anchorMatrix(1,pointer),xData0(2,pointer),xData0(3,pointer));
                    end
                end
            end
        end
        xMatrix0 = [xMatrix0(:,1:noOfSensors), anchorMatrix];
    end
% <--- placing anchors on grids
% placing anchors randomly --->
elseif anchorType == 2
    xMatrix0 = xbd*rand(sDim,noOfSensors+noOfAnchors);
% <--- placing anchors randomly
% placing anchors on the coordinate axes ---> 
elseif anchorType == 3
    noOfAnchors = sDim + 1;
    xMatrix0 = xbd*rand(sDim,noOfSensors+noOfAnchors);
    anchorMatrix = sparse(sDim,noOfAnchors); 
    pointer = 1; 
    anchorMatrix(:,pointer) = zeros(sDim,1); 
    for i=1:sDim
        pointer = pointer+1; 
        anchorMatrix(:,pointer) = zeros(sDim,1);
        anchorMatrix(i,pointer) = xbd*0.5;
    end
    xMatrix0 = [xMatrix0(:,1:noOfSensors), anchorMatrix];
% <--- placing anchors on the coordinate axes
% placing anchors in a small box near the center ---> 
elseif anchorType == 4
    noOfAnchors = sDim + 1;
    xMatrix0 = xbd*rand(sDim,noOfSensors+noOfAnchors);
    pointer = 1;
    anchorMatrix(:,pointer) = 0.5*xbd*ones(sDim,1); 
    for i=1:sDim
        pointer = pointer+1; 
        anchorMatrix(:,pointer) = 0.5*xbd*ones(sDim,1);
        anchorMatrix(i,pointer) = 0.6*xbd;
    end
    xMatrix0 = [xMatrix0(:,1:noOfSensors), anchorMatrix];
elseif anchorType == 10
    noOfAnchors = 0;    
    xMatrix0 = xbd*rand(sDim,noOfSensors);    
end
% <--- placing anchors in a small box near the center

% full(xMatrix0)

randSeed = randSeed + 13;

[distanceMatrix0] = computeDistance(degreeLB1,degreeLB2,noOfSensors,xMatrix0,radiorange,noisyFac,randSeed); 

fprintf('## elapsed time for generating a sensor network problem = %8.2f\n',toc(startingTime)); 

if 0 == 1
    sensorDistMat = distanceMatrix0(:,1:noOfSensors)+distanceMatrix0(:,1:noOfSensors)'+(noOfSensors+1)*speye(noOfSensors,noOfSensors);
    permutation = symamd(sensorDistMat);
    UMat = chol(sensorDistMat(permutation,permutation));
    figure(1);
    spy(UMat + UMat');
    % full(distanceMatrix0)
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [distanceMatrix0] = computeDistance(degreeLB1,degreeLB2,noOfSensors,xMatrix0,radiorange,noisyFac,randSeed)

% global ORIGINAL
% 
% if exist('ORIGINAL','var') ~= 1
%     fprintf('\n\n## Define ORIGINAL = 0 or 1 in the workspace.\n');   
%     fprintf('If you want to use orignal functions by Kojima sensei, type the following command in the workspace:\n');
%     fprintf('>> global ORIGINAL\n');
%     fprintf('>> ORIGINAL = 1;\n\n');
%     fprintf('If you want to use new functions modified by Waki, type the following command in the workspace:\n');
%     fprintf('>> global ORIGINAL\n');
%     fprintf('>> ORIGINAL = 0;\n\n');   
%     error('Please retry.'); 
% end
% 
% 
% if ORIGINAL == 0
%     fprintf('\n## Solve this problem by new functions.\n\n');
% elseif ORIGINAL == 1
%     fprintf('\n## Solve this problem by original functions.\n\n');
% else
%     fprintf('\n\n## Define ORIGINAL = 0 or 1 in your Command Window.\n');   
%     fprintf('If you want to use orignal functions by Kojima sensei, type the following command in your Command Window:\n');
%     fprintf('>> global ORIGINAL\n');
%     fprintf('>> ORIGINAL = 1;\n\n');
%     fprintf('If you want to use new functions modified by Waki, type the following command in your Command Window:\n');
%     fprintf('>> global ORIGINAL\n');
%     fprintf('>> ORIGINAL = 0;\n\n');   
%     error('Please retry.'); 
% end

sDim = size(xMatrix0,1);
% degreeLB0 = 1.0e10; 
% degreeBound = (sDim+1)*5; 
randn('seed',randSeed);
noOfAnchors = size(xMatrix0,2) - noOfSensors;
distanceMatrix0 = sparse(noOfSensors,noOfSensors+noOfAnchors);
%countVector = sparse(1,noOfSensors);
countVector = zeros(1,noOfSensors);
for r = noOfSensors+1:noOfSensors+noOfAnchors
    for p=1:noOfSensors
        if countVector(p) < degreeLB1
            d0 = norm(xMatrix0(:,p)-xMatrix0(:,r));
            if d0 <= radiorange
                rate = max([1+randn(1,1)*noisyFac,0.1]);
                distanceMatrix0(p,r) = d0*rate;
                countVector(p) = countVector(p)+1;
            end
        end
    end
end
%     if ORIGINAL == 1
%         for q = 2:noOfSensors
%             if countVector(q) < (sDim+1)*5
%                 p = 0;
%                 while p < q-1
%                     p = p+1;
%                     d0 = norm(xMatrix0(:,p)-xMatrix0(:,q));
%                     if d0 <= radiorange
%                         rate = max([1+randn(1,1)*noisyFac,0.1]);
%                         distanceMatrix0(p,q) = d0*rate;
%                         countVector(p) = countVector(p)+1;
%                         countVector(q) = countVector(q)+1;
%                     end
%                     if countVector(q) >= (sDim+1)*5
%                         p = q;
%                     end
%                 end
%             end
%         end
%
%     elseif ORIGINAL == 0
idx = find(countVector < degreeLB1+degreeLB2);
zeroPone = repmat(0.1,1,noOfSensors+noOfAnchors);
for q = idx
    tmpIdx = repmat(q,1,q-1);
    tmpMat = xMatrix0(:,1:q-1)-xMatrix0(:,tmpIdx);
    tmpMat = tmpMat.*tmpMat;
    tmpMat = sum(tmpMat,1);
    tmpMat = sqrt(tmpMat);
    I = find(tmpMat <= radiorange);
    if ~isempty(I)
        s = degreeLB1+degreeLB2 - countVector(q);
        if s < length(I)
            I = I(1:s);
        else
            s = length(I);
        end
        tmp = [1+randn(1,s)*noisyFac;zeroPone(1,1:s)];
        tmpMat(1,I) = tmpMat(1,I).* max(tmp,[],1);
        distanceMatrix0(I,q) = tmpMat(1,I)';
        countVector(I) = countVector(I) + 1;
        countVector(q) = countVector(q) + s;
    end
end
return
