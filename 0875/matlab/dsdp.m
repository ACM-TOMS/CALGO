%*****************************************************************************
% DSDP5:  Dual-Scaling Algorithm for Positive Semidefinite Programming
% Copyright (c) 2004 by
% S. J. Benson, Y. Ye
% Last modified: 20 August 2004
%*****************************************************************************
%
% > DSDP(b,AC) attempts to solve the positive semidefinite programs
%
%  (P)  MINIMIZE trace(C*X) SUCH THAT trace(A_i*X) = b_i, i=1,...,m and X >= 0
%  (D)  MAXIMIZE   dot(b,y) SUCH THAT  C - SUM A_i*y_i >= 0.
%
%      using a dual-scaling interior-point algorithm. 
%
%      The first argument b is a dense column vector of length m.
%
%      The second argument specifies the data A and C.  Since this data usually
%      has a block structure, this argument is a p x 3 cell array, where p
%      is the number of blocks.  Each row of the cell array specifies the 
%      type, dimension, and data of one block.  
%
%      Semidefinite Cone:
%      Let A_j1 .. A_jm and C_j be the parts of the full operators associated 
%      a semidefinite block.  This data can be inserted into row j of a DSDP
%      formatted cell array by following the example below: 
%         > ACj = [ dvec(A_j1) ... dvec(A_jm) dvec(C_j)];
%         > AC{j,1}='SDP'; AC{j,2}=size(C_j,1); AC{j,3}=ACj;
%      The first cell element identifies this block as a semidefinite block.
%      The second cell element identifies n, the dimension of the block, and
%      the third cell contains the data.  Each column in the third element
%      contains one constraint matrix A_j1, .., A_jm, or C_j in sparse vector 
%      form ( See DVEC() ).
%      The number of columns in ACj is LENGTH(b)+1, and the number of rows is
%      n*(n+1)/2 where n is the number of rows and columns in C_j.
%
%      Multiple Semidefinite Cones:
%      Each semidefinite cone can be specified on a row of the cell array.
%      However, these blocks can also be grouped together.   To group
%      these blocks together, the second cell entry must be an array
%      of integers representing the number of blocks and dimension of each
%      block.  The data from the blocks should be concatenated such that
%      the number of rows increases while the number of columns remains constant.
%      For instance of AC1 and AC2 represent the semidefinite blocks with
%      dimension 4 and 7, respectively, then
%         > AC{j,1}='SDP'; AC{j,2}=[4 7]; AC{j,3}=[AC1; AC2];
%
%      Fixed Variables:
%      Fixed variables are constraints where a variable y is directly 
%      set to a value.  They can be eliminated from the problem or 
%      identified to the solver.
%      If a set of variables [yi1 ... yil] is fixed to values [fi ... fl],
%      Pass the data into a separate row of the cell array.  i.e.
%         > AC{j,1}='FIXED'; AC{j,2}=[i1 ... il ]; AC{j,3}=[f1 .. fl];
%     
%      LP Cone:
%      Scalar inequalities, or matrix inequalities of dimension 1, can be
%      grouped together in a block.  Consider dual linear constraints of 
%      the form  A'*y <= c.  They can be specified by
%         > AC{j,1}='LP'; AC{j,2}=length(c); AC{j,3}=[A' c];
%      The string 'LP' represents scalar inequalities.  The 
%      second column in the cell array represents the number of scalar 
%      inequalities in the block, and the third element provides 
%      the data.  
%
% > DSDP(b,AC,OPTIONS) specifies some options for the solver
%      See DOPTIONS;
%
% > DSDP(b,AC,OPTIONS,y0) specifies an initial solution y0 in (D).
%
%
% > [STAT,y] = DSDP() returns a structure STAT with relevant information 
%              concerning the performance of the solver and an approximate
%              solution y .  See DSTAT for more description of statistics.
%
% > [STAT,y,X] = DSDP() returns a structure containing relevant statistics, 
%                and approximate solutions to (D) and (P). The solution
%                X is a cell array with the same number of columns as the
%                cell array AC and has the same block structure as AC.  Each 
%                element of the cell array X is a vector.  
%                For semidefinite blocks, the DMAT() operator may convert
%                the solution X into a matrix.
%
%  OUTPUT has columns:
%  Iter  PP Objective      DD Objective     PInfeas  DInfeas     Mu     StepLength   Pnrm
%  --------------------------------------------------------------------------------------
%  that represent the iteration number, objective values to (PP) and (DD),
%  infeasiblity to (P) and (D), barrier  parameter, 
%  steplength for (P) and (D), and a norm representing proximity to 
%  the target on the central path.
%    
%
%  See also: DVEC, DMAT, DSPARSE, DOPTIONS, DSTAT, DERROR
%
%  Examples: MAXCUT, GPP, THETAPROBLEM, CONTROL, ETP, READSDPA, READSEDUMI
%
% DSDP5.0 
% Copyright (c) 2004 by
% S. Benson and Y. Ye
% Last modified: August 2004
%*****************************************************************************


