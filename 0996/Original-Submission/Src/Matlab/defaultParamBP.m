function paramBP = defaultParamBP()
%%DEFAULTPARAMBP Setting default parameters for BP.m (BP Algorithm)
% Usage: 
%     paramBP = defaultParamBP();
% 
% paramBP has the following fields:
%     paramBP.maxtimeBP = 20000; % the maximum execution time for BP Algorithm
%     paramBP.maxiterBP = 40; % the maximum iteration for BP Algorithm 
%     paramBP.maxiterAPGR = 20000; % for maximul iteration for APGR.m (APGR Algorithm)
% 
%     paramBP.delta = 0; % the absolute terminal length
%     paramBP.delta1 = 1.0e-4; % the relative terminal length of the interval
%     
%     paramBP.printyes = 2; % print level \in {0,1,2,3}

    paramBP.sparseSW = 0;
    
    paramBP.maxtimeBP = 20000; % the maximum execution time for BP Algorithm
    paramBP.maxiterBP = 40;
    paramBP.maxiterAPGR = 20000;

    paramBP.delta = 0; % the absolute terminal length
    paramBP.delta1 = 1.0e-4; % the relative terminal length of the interval
    
    paramBP.printyes = 2; % print level \in {0,1,2,3}

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % the parameters below are also used, but their change might destroy 
    % the efficiency and the effectiveness of BP Algorithm!
    paramBP.delta2 = max(paramBP.delta,paramBP.delta1)/10;
    % If the relative increase of y0^k + rho t^k is less than paramBP.delta2 
    % and the relative length of the interval gets smaller than
    % paramBP.delta1, the iteration stops.

    paramBP.validSW = 1;
    paramBP.weightLBvsUB = 0.5;
    
    paramBP.scale_data = 1;
%   logical value to determine whether BP scales
%   coefficient matrices *at first*
%   for the purpose of numerical stability. (default: true)

    paramBP.Gscale_yes = 1;
%   logical value to determine whether BP scales
%   coefficient matrices *adaptively*
%   for the purpose of numerical stability. (default: true)

    paramBP.validTOL = sqrt(10)*1e-12; % epsilon
    paramBP.tol = 1e-12; % delta
    
    paramBP.gamma = 100;
    paramBP.heuristicFISTA = true;
    paramBP.stabilize = true;
    paramBP.Linit = 0.8;
    paramBP.newBPswitch = false;
    paramBP.betterInit = true;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end 