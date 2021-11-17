function [sol, info] = BBCPOP(objPoly, I01, Icomp, relaxOrder, params)
%BBCPOP computes a valid lower bound of a POP
%   Usage:
%       [sol, info] = BBCPOP(objPoly, I01, Icomp, relaxOrder, params);
%
%   This function converts 
%   the binary, box, complementarity constrained polynomial optimization problem (BBCPOP):
%
%      zeta* = min. evalPoly(objPoly,x)
%      s.t. x_i \in [0,1]   ( i \in Ibox )
%           x_j \in {0,1}   ( j \in I01 )
%           x^\alpha = 0    ( \alpha = Icomp(ii,:) )
%
%   into the following DNN (or SDP) relaxation problem:
%
%      min  Q0vec'*x
%      s.t. H0vec'*x = 1
%           trace(x) \leq UbdIX
%           x \in PSDcone  (positive semidefinite cone)
%           x \in polyCone (polyhedral cone). 
%
%   Then BP Algorithm is applied to the dual 
%
%       max  y0 + UbdIX*y1
%       s.t. Q0vec - H0vec*y0 + I*y1 = Y1 + Y2, y1 <= 0
%            Y1 \in (PSDcone)^*, Y2 \in (polyCone)^*.
% 
%   Input format:
%      objPoly: polynomial in sparsePOP format, i.e.,
%               struct which contains the following fields
%         .supports : 2-dimensional array whose 
%                     ii-th row is the support of the ii-th monomial
%                     (jj-th column corresponds to the the jj-th variable).
%         .coef     : column vector whose ii-th element is 
%                     the coefficient of the ii-th monomial.
%         (.typeCone .sizeCone .dimVar .degree .noTerms are not needed.)
%
%      I01    : logical row vector which specifies binary variables.
%      Icomp  : logical 2-dimensional array whose ii-th row is
%               the support of ii-th complementarity condition.
%      relaxOrder: relaxation order (integer value).
%      params : struct to specify options. It can have the following fields:
%           .sparseSW  : boolean that determines whether to exploit sparsity 
%                        false for dense DNN relaxation
%                        true for sparse DNN relaxation (default)
%           .maxtimeBP : the maximum execution time for BP Algorithm
%                        20000 (default)
%           .maxiterAPGR :  the maximum iteration for APGR Algorithm
%                           20000 (default)
%           .delta1    : the relative tolerance for BP Algorithm
%                        1e-4 (default)
%                        set 0 if delta2 is used
%           .delta2    : the absolute tolerance  for BP Algorithm
%                        0 (default)
%                        set delta2 in (0,1) if zeta* is integer
%                        set 0 otherwise
%
%   Output format:
%      sol:
%           .y0init     : the initial point for BP Algorthm
%           .LBv        : the valid lower bound computed
%           .LB         : the lower bound computed
%           .UB         : the upper bound computed
%           .Y1         : \in K_1, the approximate optimal solution computed
%           .Y2         : \in K_2, the approximate optimal solution computed
%     info:
%           .iterBP     : the number of iterations in BP Algorithm
%           .timeBP     : the execution time of BP Algorithm
%           .iterAPGR   : the total number of iterations in APGR Algorithm
%           .termcodeBP : the terminal code of BP Algorithm
%
% Reference:
% N. Ito, S. Kim, M. Kojima, A. Takeda, and K.-C. Toh.
% A Sparse Doubly Nonnegative Relaxation of Polynomial Optimization Problems
% with Binary, Box and Complementarity Constraints,
% Research Rport B-48?, Department of Mathematical and Computing Sciences, 
% Oh-Okayama, Meguro-ku, Tokyo 152-8552, March 2018. 

[Q0vec, H0vec, PSDcone, polyCone, UbdIX] ...
    = BBCPOPtoDNN(objPoly, I01, Icomp, relaxOrder, params);
% params.newBPswitch = true;
[sol, info] =...
    solveDNN(Q0vec, H0vec, PSDcone, polyCone, UbdIX, params);
end

