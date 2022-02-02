function [param,SDPobjValue,POP,cpuTime,SeDuMiInfo,SDPinfo] = ... 
    SDPrelaxationMex(param,objPoly,ineqPolySys,lbd,ubd)

% 
% SDPrelaxationMex wroks as SDPrelaxation to solve a polinomial optimization 
% problem described in SparsePOP format. The main part of SDPrelaxationMex 
% is written in C++ to speed up the cnversion of a given POP into an SDP 
% relaxation problem. 
% Usage: 
% [param,SDPobjValue,POP,cpuTime,SeDuMiInfo,SDPinfo] = ...
%		SDPrelaxation(param,objPoly,ineqPolySys,lbd,ubd);
%
% NOTICE: This program does not check the validness of inputs. Users must 
% check it by subPrograms/Mfiles/checkPOP.m in advance when they use this program.
% Moreover, input polynomials must not have monomials whose coefficients are zero. 
%
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
SDPobjValue = [];
cpuTime.SeDuMi = 0.0;
SeDuMiInfo = [];
SDPinfo.rowSizeA = 0;
SDPinfo.colSizeA = 0;
SDPinfo.nonzerosInA = 0;
SDPinfo.noOfFreevar = 0;
SDPinfo.noOfLPvariables = 0;
SDPinfo.SOCPblock = [];
SDPinfo.SDPblock = [];
if issparse(lbd)
	lbd = full(lbd);
end
if issparse(ubd)
	ubd = full(ubd);
end

%************************************************************************
fileId = 0;
if ischar(param.detailedInfFile)
	fileId = fopen(param.detailedInfFile,'a');  
end

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

%*************************************************************************
%                       Make special data-type used in mex function
%*************************************************************************
[   typelist,...
    sizelist,...
    degreelist,...
    dimvarlist,...
    notermslist,...
    supdata,...
    coefdata ] = make_mexdata(objPoly, ineqPolySys);

%[s,w]=system('top -b -n 1 | grep MATLAB | head -1 ');
%if s == 0
%	disp(w);
%end
%************************************************************************
%                        Convert POP (Part1)
%************************************************************************
%fprintf('+++ Convert Part1 +++\n');
[   new_typelist,...
    new_sizelist,...
    new_degreelist,...
    new_dimvarlist,...
    new_notermslist,...
    new_lbd,...
    new_ubd,...
    new_supdata,...
    new_coefdata,...
    objconst,...
    scalevalue,...
    bvect,...
    permmatrix,...
    cspmatrix,...
	lbdIdx,...
    ubdIdx] = mexconv1( typelist, sizelist, degreelist, dimvarlist,notermslist,lbd,ubd,supdata,coefdata,param);

%[s,w]=system('top -b -n 1 | grep MATLAB | head -1 ');
%if s == 0
%	disp(w);
%end
% transform objPoly into sparsePOP format.
NewobjPoly.typeCone = new_typelist(1);
NewobjPoly.sizeCone = new_sizelist(1);
NewobjPoly.degree = new_degreelist(1);
NewobjPoly.dimVar = new_dimvarlist(1);
NewobjPoly.noTerms = new_notermslist(1);
NewobjPoly.supports = sparse(new_supdata(:,1:NewobjPoly.noTerms)');
NewobjPoly.coef = new_coefdata(1,1:NewobjPoly.noTerms)';
%full([NewobjPoly.supports, NewobjPoly.coef])
if exist('mexconv4') == 3
[   typelist,...
    sizelist,...
    degreelist,...
    dimvarlist,...
    notermslist,...
    supdata,...
    coefdata ] = make_mexdata(NewobjPoly, ineqPolySys0);
end
%[s,w]=system('top -b -n 1 | grep MATLAB | head -1 ');
%if s == 0
%	disp(w);
%end


%************************************************************************
%              Chordal Extension by Cholesky decompostion
%************************************************************************
if param.sparseSW == 1
	cspmatrix = cspmatrix + cspmatrix';
	cspmatrix = cspmatrix + 5*new_dimvarlist(1)*speye(new_dimvarlist(1));
	[oriIdx,stats] = symamd(cspmatrix);
	[extmatrix,p] = chol(cspmatrix(oriIdx,oriIdx));
	if(p > 0)
		error('Correlative sparsity matrix is not positive definite\n');
	end
	extmatrix = extmatrix';
else
    s = size(cspmatrix,1);
    extmatrix = speye(s);
    oriIdx = (1:s);
end
if ~issparse(extmatrix)
	extmatrix = sparse(extmatrix);
end
if ~issparse(oriIdx)
	oriIdx = sparse(oriIdx);
end
%[s,w]=system('top -b -n 1 | grep MATLAB | head -1 ');
%if s == 0
%	disp(w);
%end


%************************************************************************
%                           Convert POP (Part2)
%************************************************************************
%fprintf('+++ Convert Part2 +++\n');
[ SDPA , linearterms] = mexconv2(...
    new_typelist, ...
    new_sizelist, ...
    new_degreelist, ...
    new_dimvarlist, ...
    new_notermslist, ...
    new_lbd,...
    new_ubd,...
    new_supdata,...
    new_coefdata,...
    param,...
    extmatrix,...
    oriIdx,...
    lbdIdx,...
	ubdIdx);

% PSDP is converted into the dual problem.  
[A, b, c, K] = SDPAtoSeDuMi(SDPA); 

% This command works well on Linux and Unix.
%[s,w]=system('top -b -n 1 | grep MATLAB | head -1 ');
%if s == 0
%	disp(w);
%end

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

%
% In the case of relax order = max(deg(f_i)), when applying the reducing
% technique, SparsePOP sometimes detects that the SDP problem obtained
% from its Sum of Squares problem is infeasible. In this case, By increasing 
% relax order or adding the valid inequalities, a user can overcome
% the ineasibility. 
%
if length(linearterms) ~= objPoly.dimVar && param.reduceMomentMatSW == 1
    fprintf('The dual SDP is infeasible.');
    cpuTime.conversion = cputime - startingTime1;
    cpuTime.total = cpuTime.conversion;
    SeDuMiInfo.dinf = 1;
    return
end

cpuTime.conversion = cputime - startingTime1;

%************************************************************************
%                          SeDuMi(solve SDP)
%************************************************************************

%SeDuMi
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
        %****************************************************************
        %Generate Solutions
        %****************************************************************
        trans.Amat = permmatrix;
        trans.bVect = bvect;
        trans.objConstant = objconst;
        trans.objValScale = scalevalue(1);
        % computing an approximate solution of the POP
       	[POP,SDPobjValue] = genApproxSolution(y,...
     			linearterms,SDPobjValue,NewobjPoly,objPoly0,...
			ineqPolySys0,lbd0,ubd0,param,trans,fixedVar);
    end
elseif (SDPinfo.infeasibleSW == -1)
	yVect = zeros(rowSizeAorg,1);  
	yVect(nzRowIdxUMat,1) = y;
	y = PMat'*yVect;
	SDPobjValue = -c'*x;
	trans.Amat = permmatrix;
	trans.bVect = bvect;
	trans.objConstant = objconst;
	trans.objValScale = scalevalue(1);
	% computing an approximate solution of the POP
       	[POP,SDPobjValue] = genApproxSolution(y,...
     	    	   	linearterms,SDPobjValue,NewobjPoly,objPoly0,...
			ineqPolySys0,lbd0,ubd0,param,trans,fixedVar);
elseif param.SeDuMiSW == 0 && fileId > 0
	writeSeDuMiInputData(fileId,0,A,b,c,K);
end
if fileId > 0
	fclose(fileId);
end
cpuTime.total = cputime -startingTime1;

%fprintf('.SparsePOP 2.00 succeeded \n');
return

