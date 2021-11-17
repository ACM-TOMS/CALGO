function [sol, info] =...
    solveCOP(Q0vec, H0vec, Qpvecs, PSDcone, polyCone, params)
%SOLVECOP solves a conic optimization problem (COP) below
%   Usage:
%      [sol, info] = SOLVECOP(Q0vec, H0vec, Qpvecs, PSDcone, polyCone, params);
%
%   COP: 
%      min  Q0vec' * x
%      s.t. H0vec' * x = 1
%           trace(x) \leq UbdIX (optional)
%           x \in PSDcone (positive semidefinite cone)
%           x \in polyCone (polyhedral cone).
% 
%   The dual of this problem is: 
%       max  y0 + UbdIX*y1
%       s.t. Q0vec - H0vec*y0 + I*y1 = Y1 + Y2, y1 <= 0
%            Y1 \in (PSDcone)^*, Y2 \in (polyCone)^*.
%   Input:
% 
%      Q0vec, H0vec, PSDcone, polyCone:
%        The output of BBPOPtoCOP. See 
%        >> help BBCOPtoCOP
%        for details.
%      Qpvecs is used to solve an equivalent COP by 'sdpnalplus', 'sedumi',
%      'sdpt3', 'sdpa' and 'mosek'
%
%      params: struct to specify options. It can have the following fields:
%          .solver: string value to determine the solver/algorithm. (default: 'BP').
%                   Other possible values are 'sdpnal', 'sdpt3', 'sedumi', 'sdpa'.
%                   The corresponding software is required to be installed.
%          
%          .UbdIX : An upper bound of the sum of the trace of the moment matrices.
%                   It is used for computing a valid lower bound.
%                   ( default: sum(PSDcone) )
%
%          The following option is only valid for 'BP' and 'sdpnal'.
%          .maxtime        : The maximum computation time. (default 20000)
%          
%          The following options are only valid for 'BP'
%          .maxiterBP      : The maximum number of iteration of Bisection. (default: 40)
%          .delta          : tolerance for bisection (default: 0)
%          .delta1         : relative tolerance for bisection (default: 1e-6)
%          .maxiterAPGR    : The maximum number of iteration of FISTA. (default: 20000)
%          
%          .scale_data     : logical value to determine whether BP scales
%                            coefficient matrices *at first*
%                            for the purpose of numerical stability. (default: true)
%          .Gscale_yes     : logical value to determine whether BP scales
%                            coefficient matrices *adaptively*
%                            for the purpose of numerical stability. (default: true)
%          .validTOL       : tolerance for feasibility. (default: 3e-12)
%          .tol            : tolerance for optimality. (default: 1e-12)
%          .heuristicFISTA : logical value to determine whether FISTA uses
%                            heuristic stopping criteria. (default: true)
%                            Setting it to true may deteriorates the accuracy of BP. 


    params = fillUnspecifiedParams(params, defaultParamsSolver());
    info = [];
    
    if ~isfield(params, 'UbdObjVal')
        params.UbdObjVal = sum(Q0vec(Q0vec>0));
    end
    if ~isfield(params, 'UbdIX')
        params.UbdIX = sum(polyCone.varStructure.sizeblk);
    end

    K.s = PSDcone;
    projMap.K1     = @(x) blkprojSDPBP(x, K);
    projMap.K1Star = @(x) blkprojSDPBP(x, K);
    projMap.K2     = @(x) projOntoPolyConeMex(x, polyCone);
    projMap.K2Star = @(x) x + projOntoPolyConeMex(-x, polyCone);

    if strcmp(params.solver, 'BP')
        UB = params.UbdObjVal;
        LB = full(sum(Q0vec(Q0vec<0)));
        [sol, infoBP] =...
            BP(Q0vec, H0vec, K, projMap, LB, UB, params);
        info.iterBP = infoBP.iter; 
        info.timeBP = infoBP.timeBisection; 
        info.iterAPGR = infoBP.APGiter;
        info.termcodeBP = infoBP.break_yes;
        return
    end
    
    if strcmp(params.solver, 'sdpnal') || strcmp(params.solver, 'sdpnalplus')
        [x,y,info] = COPsdpnalplus(Q0vec, H0vec, Qpvecs, PSDcone, polyCone, params);
        sol.x = x;
        sol.y = y;
        sol.LB = info.objval(1);
        sol.UB = info.objval(2);
        if isfield(info, 'LBv')
            sol.LBv = info.LBv;
        else
            sol.LBv = sol.LB;
        end
        return
    end

    [A,b,c,K] = COPtoSedumi(Q0vec, H0vec, Qpvecs, PSDcone, polyCone, params);
    if strcmp(params.solver, 'sedumi')
        [x,y,infoSedumi] = sedumi(A,b,c,K);
        info.cputime = infoSedumi.cpusec;
        info.iter = infoSedumi.iter;
        info.termcode = infoSedumi.numerr; 
    elseif strcmp(params.solver, 'sdpt3')
        [x,y,info] = sdpt3Sedumi(A,b,c,K);
    elseif strcmp(params.solver, 'sdpa')
        [x,y,info] = sdpaSedumi(A,b,c,K);
    elseif strcmp(params.solver, 'mosek')
        [x,y,info] = mosekSedumi(A,b,c,K);
    end
    sol.x = x;
    sol.y = y;
    sol.sedumiA = A;
    sol.sedumib = b;
    sol.sedumic = c;
    sol.sedumiK = K;
    if strcmp(params.COPrep,'basis')
        sol.LB = -x'*c;
        sol.UB = sol.LB;
    else
        sol.LB = x'*c;
        sol.UB = sol.LB;
    end
    sol.LBv = []; % sol.LB;
    
end

