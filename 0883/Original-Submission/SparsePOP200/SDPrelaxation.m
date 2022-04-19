function [param,SDPobjValue,POP,cpuTime,SeDuMiInfo,SDPinfo] = ...
		SDPrelaxation(param,objPoly,ineqPolySys,lbd,ubd)
% 
% SDPrelaxation
% solves a polinomial optimization problem described in SparsePOP format.
%
% Usage: 
% [param,SDPobjValue,POP,cpuTime,SeDuMiInfo,SDPinfo] = ...
%		SDPrelaxation(param,objPoly,ineqPolySys,lbd,ubd);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% objPoly, inEqPolySys, lbd, ubd form the SparsePOP format.
% param is a set of parameters. See below for the details.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% param: parameters actually used.
% SDPobjValue: the optimal value of SDP
% POP: 
%   POP.xVect: an approximate solution for POP
%   POP.objValue: the objective value of xVect
%   POP.absError: an absolute error
%   POP.scaledError: an scaled error
% cpuTime:
%   cpuTime.conversion: cpu time consumed to convert POP to SDP relax. 
%   cpuTime.SeDuMi: cpu time consumed by SeDuMi.
%   cpuTime.Total: total cpu time.
% SeDuMiInfo:
%   SeDuMiInfo.numerr: info.numerr of SeDuMi
%   SeDuMiInfo.pinf: info.pinf of SeDuMi
%   SeDuMiInfo.dinf: info.dinf of SeDuMi
% SDPinfo:
%   SDPinfo.rowSizeA: the number of rows of A.
%   SDPinfo.colSizeA: the number of columns of A.
%   SDPinfo.nonzeroInA: the number of nonzero elements in A.
%   SDPinfo.noOfLPvariables: the number of LP variables in SDP.
%   SDPinfo.noOfFRvariables: the number of Free variables in SDP.
%   SDPinfo.SDPblock: the row vector of sizes of SDP blocks.
%
% See UserGuide.pdf and/or the following reference:
%
% H. Waki, S. Kim, M. Kojima and M. Muramatsu,
% "Sums of Squares and Semidefinite Programming Relaxations 
% for  Polynomial Optimization Problems with Structured Sparsity",
% SIAM Journal on Optimization Vol.17 (1) 218-242 (2006).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is a component of SparsePOP 
% Copyright (C) 2007 SparsePOP Project
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

startingTime1 = cputime;

%% Outputs
SDPobjValue = [];
POP.xVect = [];
POP.objValue = [];
POP.absError = [];
POP.scaledError = [];
cpuTime.SeDuMi = 0.0;
SeDuMiInfo = [];
SDPinfo.rowSizeA = 0;
SDPinfo.colSizeA = 0;
SDPinfo.nonzerosInA = 0;
SDPinfo.noOfFreevar = 0;
SDPinfo.noOfLPvariables = 0;
SDPinfo.SOCPblock = [];
SDPinfo.SDPblock = [];

% saving the original POP information
[objPoly0,ineqPolySys0,lbd0,ubd0] = saveOriginalPOP(objPoly,...
    ineqPolySys,lbd,ubd); 

% Kojima 06/11/2007 --->
%
% param.reduceAMatSW = 1 ---> Reducing the coeeficient matrix A of SeDuMi;
%                    = 0 ---> Not reducing the coeeficient matrix A of SeDuMi;
%
infeasibleSW = 0;
% infeasibleSW = 2 ---> POP is found infeasible in nDimZero2 before applying SeDuMi; 
% infeasibleSW = 1 ---> Primal SDP is found infeasible in solveSquareSDP before applying SeDuMi;  
% infeasibleSW = 0 ---> apply the sparse SDP relaxation; 
% infeasibleSW = -1 ---> SDP is solved by solveSquareSDP and has a unique
%                        feasible solution; SeDuMi is not applied; 
% infeasibleSW = -2 ---> POP is found to have a unique feasible solution nDimZero2;
%
% SDPinfo.infeasibleSW = infeasibleSW, later; 
%
% SDPinfo.reduceAMatSW = 2; all the variables of POP are fixed in deleteVar; 
% SDPinfo.reduceAMatSW = 1; some equality constraints of the SeDuMi format primal SDP are linearly dependent; 
% SDPinfo.reduceAMatSW = 0; no variable of POP is fixed; 
%
% Removing redundant variables
deleteVarSW = 1;
if (param.reduceAMatSW == 1) && (deleteVarSW == 1)
    [objPoly,ineqPolySys,lbd,ubd,fixedVar] = deleteVar(objPoly,ineqPolySys,...
        lbd,ubd,param);
%	for i=1:size(ineqPolySys,2)
%		full([ineqPolySys{i}.supports, ineqPolySys{i}.coef])
%	end
else
    fixedVar = [];
end

% If all variables vanish from the original POP by "deleteVar.m", 
% SDP relaxations are not applied into the POP.
if objPoly.dimVar == 0
    SDPinfo.reduceAMatSW = 2; 
    % SDPinfo.reduceAMatSW = 2 ---> all the variables of POP are fixed in
    % deleteVar; 
    fprintf('## All variables of POP are fixed.\n');
    [POP, SDPobjValue,cpuTime,SeDuMiInfo,infeasibleSW] = nDimZero2(objPoly,objPoly0,...
        ineqPolySys0,lbd0,ubd0,POP,SDPobjValue,cpuTime,startingTime1,...
        SeDuMiInfo, fixedVar);
    SDPinfo.infeasibleSW = infeasibleSW; 
    % SDPinfo.infeasibleSW = 2  ---> POP is found infeasible in nDimZero2 
    %                                before applying SeDuMi; 
    % SDPinfo.infeasibleSW = -2 ---> POP is found to have a unique feasible 
    %                                solution in nDimZero2;
    return
elseif  ~isempty(fixedVar)
    %noFixedVar = length(fixedVar);
    noFixedVar = size(fixedVar,1);
    if noFixedVar == 1
        fprintf('## 1 variable of POP is fixed.\n');
    else
        fprintf('## %d variables of POP are fixed.\n',noFixedVar);
    end
end
% <--- Kojima 06/11/2007 

% printing information on the Polynomial SDP to be relaxed
fileId = 0; 
if ischar(param.detailedInfFile)
	fileId = fopen(param.detailedInfFile,'a+');  
end

if fileId > 0 
	writeParameters(fileId,param); 
	fprintf(fileId,'# POP to be solved.\n'); 
	writePOP(fileId,objPoly,ineqPolySys,lbd,ubd); 
end

% If x_i has a finite lower bound l_i,
% then a new variable y_i = x_i - l_i; 
% hence y_i becomes a nonnegative variable. 
trans.Amat = speye(objPoly.dimVar);
trans.bVect = sparse(objPoly.dimVar,1);
if param.symbolicMath == 1 && param.scalingSW == 1
		[objPoly,ineqPolySys,lbd,ubd,trans] = convert2(objPoly,...
            ineqPolySys,lbd,ubd);
end


% Incorprate lower and upper bounds into ineqPolySys 
% to strengthen the relaxation
[ineqPolySys,lbdIdx,ubdIdx] = boundToIneqPolySys(ineqPolySys,lbd,ubd);

% removing a constant term from objPoly
CTermIndex = find(any(objPoly.supports,2) == 0);
NTermIndex = find(any(objPoly.supports,2) ~= 0);
trans.objConstant = 0.0;  
if ~isempty(CTermIndex)
    trans.objConstant = sum(objPoly.coef(CTermIndex,1),1);
    objPoly.supports = objPoly.supports(NTermIndex,:);
	objPoly.coef = objPoly.coef(NTermIndex,:);
	objPoly.noTerms = length(NTermIndex);
end


% Conversion of the POP to be solved into a POP whose SDP relaxation
% is numerically stable.
trans.objValScale = 1.0;
variableScale = ones(1,objPoly.dimVar);
ineqValScale = ones(1,size(ineqPolySys,2));
writeTofile = 0;
% Scaling of lbd, ubd, objPoly and ineqPolySys.
% Scaling information, which we need to scale them back later,
% is available in objValScale,ineqValScale and variableScale.
if param.scalingSW == 1
    [objPoly,ineqPolySys,lbd,ubd,trans.objValScale,ineqValScale,...
        variableScale] = scalingPOP(objPoly,ineqPolySys,lbd,ubd);
    writeTofile = 1;
end
trans.Amat = trans.Amat * diag(variableScale);
if abs(param.perturbation) > 1.0e-12
    randSeed = 117;
    objPoly  = perturbObjPoly(objPoly, param.perturbation,randSeed);
    writeTofile = 1;
end
if param.eqTolerance > 1.0e-12
    ineqPolySys = relax1EqTo2Ineqs(objPoly,ineqPolySys,param.eqTolerance);
    writeTofile = 1;
end

if fileId > 0 && writeTofile == 1
    fprintf(fileId,'# Scaled and modified POP to be solved\n');
    writePOP(fileId,objPoly,ineqPolySys,lbd,ubd);
end

%{
full([objPoly.supports, objPoly.coef])
fprintf('objPoly\n');
for j=1:objPoly.noTerms
	fprintf('coef = %20.18f\n',objPoly.coef(j));
end
for i=1:size(ineqPolySys,2)
	fprintf('ineqPolySysi -- %d\n',i);
	for j=1:ineqPolySys{i}.noTerms
		fprintf('coef = %20.18f\n',ineqPolySys{i}.coef(j));
	end
end
%}

% Analyzing the correlation sparsity of POP and generating the maximal
% cliques of the csp graph induced from POP if param.sparseSW == 1. 
clique = genClique(objPoly,ineqPolySys,param.sparseSW);
if fileId > 0 
	writeClique(fileId,clique);
end

% Construction of basisIndices used in basisSupports.  
% basisIndices vary depending on param.sparseSW and
% param.multiCliquesFactor
basisIndices= genBasisIndices(objPoly,ineqPolySys,clique.Set,param);
if fileId > 0 
	writeBasisIndices(fileId,basisIndices);
end


% Construction of basisSupports. 
basisSupports = genBasisSupports(objPoly,ineqPolySys,param,basisIndices); 
if fileId > 0 
	writeBasisSupports(fileId,basisSupports);
end

[basisSupports,ineqBasis] = reduceSupSets(objPoly,ineqPolySys,...
    basisSupports,param);
[CompSup,ConstraintInfo] = separateSpecialMonomial(ineqPolySys,param);
[basisSupports,momentSup,ineqBasis] = substituteEq(basisSupports,...
    ineqBasis,ineqPolySys,CompSup,param);

if fileId > 0 
	fprintf(fileId,'# basisSupports after reduction\n'); 
	writeBasisSupports(fileId,basisSupports);
end

% Add bounds to all monoials
[ineqPolySys,basisSupports,boundList] = addBoundToPOP(ineqPolySys,...
        basisSupports,lbd,ubd,momentSup,ineqBasis,lbdIdx,ubdIdx,param);

% relaxation of the Polynomial SDP into an SDP in the SDPA format
[SDPA,xIdxVec] = PSDPtoLSDP(objPoly,ineqPolySys,basisSupports,boundList,...
    CompSup,ConstraintInfo);

%%
%% In the case of relax order = max(deg(f_i)), when applying the reducing
%% technique, we often can find the infeasibility of SDP problem obtained
%% from a Sum of Squares problem. By increasing relax order or adding the
%% valid inequalities, a user can overcome the infeasibility.
%%
Degree = sum(xIdxVec,2);
linearterms = find(Degree==1);
linearterms = linearterms -1;
if length(linearterms) ~= objPoly.dimVar && param.reduceMomentMatSW == 1
    fprintf('The dual SDP is infeasible.');
    SeDuMiInfo.dinf = 1;
    cpuTime.conversion = cputime - startingTime1;
    cpuTime.total = cpuTime.conversion;
    return
end

if isfield(param,'sdpaDataFile') && ~isempty(param.sdpaDataFile)
    sdpaDataId = fopen(param.sdpaDataFile,'w');
    writeSDPAformatData(sdpaDataId,param,SDPA);
    fclose(sdpaDataId);
end

% PSDP is converted into the dual problem.

[A, b, c, K] = SDPAtoSeDuMi(SDPA);

% Kojima 06/11/2007 --->
% Cheking the linear dependence of the row vectors of the matrix A.
% If they are linearly dependent, then  
% (a) reduce the system of equations A x = b,
% (b) detect whether it is infeasible, and/or
% (c) compute the unique feasible solution if the system is nonsingular and
% square.
% Set reducedMatSW = 0 to suppress this function
%
%
if param.reduceAMatSW == 1
    % infeasibleSW =  1 ---> primal SDP is infeasible
    % infeasibleSW =  0 ---> feasible and to be solved by SeDuMi if param.SeDuMiSW == 1
    % infeasibleSW = -1 ---> feasible, the SDP is square and to be solved by solveSquareSDP 
    %                        Not necessary to apply SeDuMi
    [rowSizeAorg,colSizeAorg] = size(A);  
    reduceAMatSW = 0;
    [reduceAMatSW,infeasibleSW,A,b,PMat,nzRowIdxUMat,x,y] = reduceSizeA4(A,b,c,K); 
%    infeasibleSW
else 
    reduceAMatSW = 0; 
end
% <--- Kojima 06/11/2007

% information on SDP solved
SDPinfo = getSDPinfo(A,K);

% Kojima 06/11/2007 --->
SDPinfo.reduceAMatSW = reduceAMatSW; 
SDPinfo.infeasibleSW = infeasibleSW;  
% <--- Kojima 06/11/2007 

cpuTime.conversion = cputime - startingTime1;

%************************************************************************
%                          SeDuMi(solve SDP)
%************************************************************************
% Kojima 06/11/2007 --->
if (param.SeDuMiSW == 1) && (SDPinfo.infeasibleSW == 0)
    [x,y,SeDuMiInfo] = solveBySeDuMi(fileId,A,b,c,K,param);
    if SDPinfo.reduceAMatSW == 1
        yVect = zeros(rowSizeAorg,1);  
        yVect(nzRowIdxUMat,1) = y;
        y = PMat'*yVect;
    end
% <--- Kojima 06/11/2007    
    cpuTime.SeDuMi = SeDuMiInfo.cpusec;
    if SeDuMiInfo.pinf == 0 && SeDuMiInfo.dinf == 0
        SDPobjValue = -c'*x;
        % computing an approximate solution of the POP
        [POP,SDPobjValue] = genApproxSolution(y,linearterms,SDPobjValue,...
            objPoly,objPoly0,ineqPolySys0,lbd0,ubd0,param,trans,fixedVar); 
    end
elseif (SDPinfo.infeasibleSW == -1)
	yVect = zeros(rowSizeAorg,1);  
	yVect(nzRowIdxUMat,1) = y;
	y = PMat'*yVect;
    	SDPobjValue = -c'*x;
	% computing an approximate solution of the POP
    [POP,SDPobjValue] = genApproxSolution(y,linearterms,SDPobjValue,...
        objPoly,objPoly0,ineqPolySys0,lbd0,ubd0,param,trans,fixedVar); 
elseif (param.SeDuMiSW == 0) && fileId > 0
	writeSeDuMiInputData(fileId,0,A,b,c,K);
end
cpuTime.total = cputime -startingTime1;
if fileId > 0
    fclose('all');
end

return

% $Header: /home/waki9/CVS_DB/SparsePOPdev/subPrograms/sparsePOPmain.m,v 1.4 2007/01/15 10:41:26 waki9 Exp $
