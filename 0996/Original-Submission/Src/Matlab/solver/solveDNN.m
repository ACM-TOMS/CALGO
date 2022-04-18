function [sol, info] = solveDNN(Q0vec, H0vec, PSDcone, polyCone, UbdIX, params)
%SOLVEDNN solves a DNN relaxation problem
%   Usage:
%      [sol, info] = SOLVEDNN(Q0vec, H0vec, Qpvecs, PSDcone, polyCone, params);
%
%      min  Q0vec' * x
%      s.t. H0vec' * x = 1
%           trace(x) \leq UbdIX
%           x \in PSDcone (positive semidefinite cone)
%           x \in polyCone (polyhedral cone). 
% 
% The dual of this problem 
%       max  y0 + UbdIX*y1
%       s.t. Q0vec - H0vec*y0 + I*y1 = Y1 + Y2, y1 <= 0
%            Y1 \in (PSDcone)^*, Y2 \in (polyCone)^*.
% is solved. 
%
%   Input:
%      Q0vec, H0vec, PSDcone, polyCone, UbdIX:
%        The output of BBCPOPtoDNN. See 
%        >> help BBCPOPtoDNN
%        for details.
%
%      params: struct to specify options. It can have the following fields:
%          .maxiterBP      : The maximum number of iteration of Bisection. (default: 40)
%          .delta          : tolerance for bisection (default: 0)
%          .delta1         : relative tolerance for bisection (default: 1e-6)
%          .maxiterAPGR    : The maximum number of iteration of FISTA. (default: 20000)
%          .maxtime        : The maximum computation time. (default)
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

params.UbdIX = UbdIX;
[sol, info] =...
    solveCOP(Q0vec, H0vec, [], PSDcone, polyCone, params);

end