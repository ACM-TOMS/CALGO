function [Q0vec, H0vec, PSDcone, polyCone, UbdIX] ...
    = BBCPOPtoDNN(objPoly, I01, Icomp, relaxOrder, params)
%%BBCPOPTODNN converts binary and box constrained polynomial optimization problem 
% into cone relaxation problem
%   
%   Usage:
%       [Q0vec, H0vec, PSDcone, polyCone] = ...
%                  BBCPOPTODNN(objPoly, I01, Icomp, relaxOrder, params);
%
%   This function converts 
%   the binary and box constrained polynomial optimization problem (BBPOP):
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
%  The dual of this problem is written as 
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
%          .sparseSW  : boolean that determines whether to exploit
%                       sparsity. (default: true)
%
%   Output format:
%      Q0vec   : vectorized coefficient matrix (sedumi format)
%      H0vec   : vectorized coefficient matrix (sedumi format)
%      PSDcone : size of each matrix variable (corresponds to K.s of sedumi format)
%      polyCone: struct which contains the following fields
%         .varStructure : struct which contains the following fields
%             .supports : supports of variables which appear in moment matrix
%             .objic    : support of ii-th term of poly is supports(objic(ii),:)
%             .momentic, .blkIdx, .rowIdx, .colIdx :
%                   support of (blkIdx(ii),rowIdx(ii),colIdx(ii)) element 
%                   is supports(momentic(ii),:). (only lowtri element)
%             .BMat     : each column is vectorized basis matrix of moment matrix
%             .sizeblk  : size of each variable matrix (equals to PSDcone)
%         .chain    : logical 2-dimensional array whose column represents a chain 
%         .nonneg   : flag of nonnegative constraints. % for SDPsolver
%         .eq0Var   : logical vector 
%         .lbd      : lower bounds of SDP variables
%         .ubd      : upper bounds of SDP variables
%     UbdIX   : the upper bound for trace(x) in the primal DNN relaxaton
%               problem above.

params.relaxOrder = relaxOrder;
[Q0vec, H0vec, ~, PSDcone, polyCone, UbdIX] ...
    = BBPOPtoCOP(objPoly, [], I01, Icomp, params);
end

