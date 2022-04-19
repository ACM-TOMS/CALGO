%
% This function relaxes a polynomial SDP into a standard SDP.
% We assume the polynomilal SDP has been generated from a polynomial
% optimization.
%
% Problem of the form
%	minimize 	f_0(\x)
%	subject to 	f_i(\x) \succeq \O (j=1,2,\ldots,m).
%
% Input --- a polynomial SDP of the following form
%	minimize 	f_0(\x)
%	subject to 	f_i(\x) \u_{\GC_j}(\x) \u_{\GC_j}(\x)^T \succeq \O (j=1,2,\ldots,m),
%				\u_{\GC_j}(\x) \u_{\GC_j}(\x)^T \succeq \O (j=m+1,\ldots,m')
%
% Here, f_0(\x) is 'objPoly', f_i(\x) is 'inEqPolySys{i}' and
% \u_{\GC_j}(\x) is 'basisSupports{j}'
%
% This function returns SDPA sparse format of relaxed SDP derived from
% given POP and all monomials in the SDP.
%
% SDPA is the structure with the following fields:
%
% SDPA.blockStruct: the block structre vector
% SDPA.nBlock:      the number of block
% SDPA.mDim:        the nubmer of the primal variables
% SDPA.bVect:       in Primal SDP, coefficients of objective function.
% SDPA.typeCone:    indicate what the block derived from
%                   real-valued nonnegative constraint,
%                   equation or SDP constraint.
% SDPA.sparseMatrix:this is matrix with 5 columns
%                   First column is the number of matrix.
%                   Second column is the number of block.
%                   Third column is the row index.
%                   Fourth column is the column index.
%                   Last column is the coefficient.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% M. Kojima, 02/06/2005
% The main function:
% 	PSDPtoLSDP(objPoly,inEqPolySys,basisSupports,ConstraintInfo,param);
% This module also inclues the following functions.
%	AddSupport(basisSupports,inEqPolySys,Row,Row);
%	makeMoment(basisSupports);
%	make_varList_basisSupports(basisSupports,nBlock,removeSup);
%	make_varList_tC1(basisSupports,inEqPolySys,nBlock,removeSup);
%	make_varList_tC3(basisSupports{j},inEqPolySys{j},nBlock,removeSup);
%	make_varList_eq(basisSupports,inEqPolySys,nBlock,removeSup);
%	removeComp(Sup,Row,Col,Coef,oldsize,removeSup);
%	SOCPtoSDP(anIneqPoly);
%	SupTerm1(inEqPolySys,ConstraintInfo,param);
%	varListtoSDPA(varList, SDPA);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [SDPA,xIdxVec] = PSDPtoLSDP(objPoly,inEqPolySys,basisSupports,boundList,CompSup,ConstraintInfo)

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

noOfInequalities = size(inEqPolySys,2);
noOfBasisSupports = size(basisSupports,2);
CompConst = ConstraintInfo.Term1.comp;
inEqIdx = setdiff((1:noOfInequalities),CompConst);

%
% varLists has 'nDim+4' columns
% varList consists of 4 pars.
%
% Frist part is varList(:,1:nDim), which has all monomials in
% Polynomial SDP.
%
% Second part varList(:,nDim+1),  is block number that the monomial
% appear in.
%
%
% Third part varList(:,nDim+2:nDim+3), is row and column indeces of
% the monomial in the block number.
%
% Last part is coefficient.
%
tempList = cell(1+noOfBasisSupports,1);
tempList{1} = [objPoly.supports,sparse(objPoly.noTerms,3),objPoly.coef]';
bSt = cell(1+noOfBasisSupports,1);
bSt{1} = [];
Atype = cell(1+noOfBasisSupports,1);
Atype{1} = [];


%%
%% Gather SDPA information of Moment matrix constructed by
%% basisSupports
%%
nBlock = 0;
pointer = 1;
if noOfBasisSupports > noOfInequalities
    for k=(noOfInequalities+1):noOfBasisSupports
        pointer = pointer + 1;
        [temp,blockSt] = ...
            make_varList_basisSupports(basisSupports{k},nBlock,CompSup);
        if ~isempty(temp)
            nBlock = nBlock + 1;
            if blockSt > -1
                Atype{pointer} = 3;
            elseif blockSt == -1
                Atype{pointer} = 1;
            end
            bSt{pointer} = blockSt;
            tempList{pointer} = temp';
        end
    end
end
%%
%% Gather SDPA information of Localize matrix constructed by
%% basisSupports and constraints
%%

for j=inEqIdx
    typeCone = inEqPolySys{j}.typeCone;
    pointer = pointer + 1;
    if typeCone == 2
        [inEqPolySys{j}.coef,inEqPolySys{j}.typeCone] = SOCPtoSDP(inEqPolySys{j});
        %if size(basisSupports{j},1) == 1 
        %    typeCone = 2;
        %else
            typeCone = 3;
        %end
    end
    %%
    %% nonnegative constraint case
    %%
    if typeCone == 1
        [temp,blockSt] =  make_varList_tC1(basisSupports{j}, ...
            inEqPolySys{j},nBlock,CompSup);
        %%
        %% SDP constraint case
        %% SOCP case (SOCP constraint in Polynomial SDP)
        %%
    elseif typeCone == 2 || typeCone == 3
        [temp,blockSt] ...
            = make_varList_tC3(basisSupports{j},inEqPolySys{j},nBlock,CompSup);
        %%
        %% equation case
        %%
    elseif typeCone == -1
        [temp,blockSt] = make_varList_eq(basisSupports{j}, ...
            inEqPolySys{j},nBlock,CompSup);
    end
    if ~isempty(temp)
        if typeCone == 1 || typeCone == -1
            typeIdx = find(blockSt);
            blockSt = blockSt(typeIdx);
            typeIdx = zeros(1,length(typeIdx));
            SDPidx = find(blockSt > 0);
            LPidx = find(blockSt == -1);
            EQidx = find(blockSt < -1);
            typeIdx(1,SDPidx) = 3;
            typeIdx(1,LPidx) = 1;
            typeIdx(1,EQidx) = -1;
            nBlock = nBlock + length(typeIdx);
            Atype{pointer} = typeIdx;
        elseif typeCone == 2
            nBlock = nBlock + 1;
            Atype{pointer} = 2;   
        elseif typeCone == 3
            nBlock = nBlock + 1;
            Atype{pointer} = 3;          
        end
        bSt{pointer} = blockSt;
        tempList{pointer} = temp';
    end
end
if ~isempty(boundList)
    nBlock = nBlock+1;
    boundList(:,end-3) = nBlock;
    pointer = pointer + 1;
    tempList{pointer} = boundList';
    bSt{pointer} = -boundList(end,end-2);
    Atype{pointer} = 1;
end

varList = [tempList{1:pointer}]';
SDPA.blockStruct = [bSt{1:pointer}];
SDPA.typeCone = [Atype{1:pointer}];
SDPA.nBlock = length(SDPA.blockStruct);
%%
%% Convert varList into SDPA sparse format.
%%
[SDPA,xIdxVec]  = varListtoSDPA(varList, SDPA);
SDPA = aggregateBlocks(SDPA);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tempList,blockStruct] = ...
    make_varList_basisSupports(basisSupports,nBlock,removeSup)

%
% Check the input.
%
if isempty(basisSupports)
    tempList = [];
    blockStruct = [];
    return;
end

%
% Construct Moment matrix.
%
m1 = size(basisSupports,1);
[MomSup, Row, Col] = makeMoment(basisSupports);
Coef = ones(length(Row),1);

%
% Remove complementarity term.
%
[NewSup,NewRow,NewCol,NewCoef,Newm1] = removeComp(MomSup,Row, ...
    Col,Coef,m1,removeSup);

%
% Output blockStruct and tempList
%
if isempty(NewSup)
    blockStruct = [];
    tempList = [];
else
    if Newm1 == 1
        %% 2005 03 06
        %% Add the part that remove u(x) = 1 case.
        %%
        if ~any(NewSup,2)
            tempList = [];
            blockStruct = [];
            return;
        end
        blockStruct = -1;
    else
        blockStruct = Newm1;
    end

    tempList = [NewSup,repmat(nBlock+1,length(NewRow),1),NewRow,NewCol,NewCoef];
end
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tempList,blockStruct] ...
    = make_varList_eq(basisSupports,inEqPolySys, nBlock,removeSup)

% This function returns information of Polynomial equations,
% that is, monomials in this Polynomial equations, their block No,
% row and column indices, and coefficients.
%
%
% <Return Value>
% 'tempList' has nDim+4 columns and 'blockstruct'
%
% tempList(:,1:nDim) is monomial that appears in the Polynomial
% equations
% tempList(:,nDim+1) is the block number that the monomial belongs.
% tempList(:,nDim+2:nDim+3) are row and column indices in the
% matrix of block number
% tempList(:,end) is for the coefficients.
%
% blockStruct is the size of Polynomial equations.
%
%

if isempty(basisSupports)
    tempList = [];
    blockStruct = [];
    return;
end
m0 = size(inEqPolySys.supports,1);
m1 = size(basisSupports,1);
sizeCone = inEqPolySys.sizeCone;

Row = repmat((1:m1)',m0,1);
%%
%% Only add basisSupport and Support of constraints because
%% basisSupport already prepared in the case of equation.
%%
Sup = inEqPolySys.supports;
tempList = [];
blockStruct = [];
for i=1:sizeCone
	Coef = inEqPolySys.coef(:,i);
	[NewSup,NewRow,NewCol,Coef] = AddSupport(basisSupports,Sup,Row,Row,Coef);
	%%
	%% Remove the term that becomes zero by substituting complementarity.
	%%
	[NewSup,NewRow,NewCol,NewCoef,Newm1] = removeComp(NewSup,NewRow,NewCol, ...
    	Coef,m1,removeSup);

	%%
	%% Arrange Output information and Remove the part whose
	%% coefficient is zero
	%%
	lenOneTerm = length(NewRow);
	blockNo = repmat(nBlock+(1:sizeCone),lenOneTerm, 1);
	blockNo = blockNo(:);
	tmpidx = repmat((1:lenOneTerm)',sizeCone,1);
	[nonZeroIdx, temp,val] = find(NewCoef(:));
	tmpidx = tmpidx(nonZeroIdx);
	blockNo = blockNo(nonZeroIdx);
	nonZeroVal = NewCoef(nonZeroIdx);
	SDPAversion = 6.2;

	if SDPAversion <= 6.2
    	if ~isempty(NewRow) || any(any(NewSup,2))
        	blockStruct = [blockStruct,-2*Newm1];
        	AddRow = NewRow + Newm1;
        	tempList = [tempList;NewSup(tmpidx,:),blockNo,NewRow(tmpidx),NewRow(tmpidx),nonZeroVal;
            	NewSup(tmpidx,:),blockNo,AddRow(tmpidx),AddRow(tmpidx),-nonZeroVal];
    	else
        	%tempList = [];
        	%blockStruct = [];
        	%return;
    	end
	else
    	%%%For new SDPA sparse format
    	%%%But, New SDPA version has not been released yet.
    	blockStruct = -Newm1;
    	tempList = [NewSup(tmpidx,:),blockNo,Row(tmpidx),nonZeroVal];
	end
end
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tempList,blockStruct] ...
    = make_varList_tC1(basisSupports,inEqPolySys,nBlock,removeSup)
%
% This function returns information of Polynomial inequalities.
% That is, monomials in this Polynomial inequalities, their block number,
% row and column indices, and coefficients.
%
%
% <Return Value>
% 'tempList' has nDim+4 columns and 'blockstruct'
%
% tempList(:,1:nDim) is monomial that appear in the Polynomial
% inequalities
% tempList(:,nDim+1) is the block No that the monomial belongs.
% tempList(:,nDim+2:nDim+3) are row and column indices in the
% matrix of block No
% tempList(:,end) is coefficients.
%
% blockStruct is the size of Polynomial inequalities.
%
%
if isempty(basisSupports)
    tempList = [];
    blockStruct = [];
    return;
end
m1 = size(basisSupports,1);
sizeCone = inEqPolySys.sizeCone;
%%
%% Construct Moment matrix from basisSupports
%%
[MomSup, Row, Col] = makeMoment(basisSupports);
%%
%% Add Moment matrix and Support of constriant
%%
Sup = inEqPolySys.supports;
tempList = [];
blockStruct = [];
for i=1:sizeCone
	Coef = inEqPolySys.coef(:,i);
	[NewSup,NewRow,NewCol,Coef] = AddSupport(MomSup,Sup,Row,Col,Coef);
	%%
	%% Remove the term that becomes zero by substituting complementarity.
	%%
	[NewSup,NewRow,NewCol,NewCoef,Newm1] = removeComp(NewSup,NewRow, ...
    	NewCol,Coef,m1,removeSup);
	%%
	%% Arrange Output information and Remove the part whose
	%% coefficient is zero
	%%
	%% Second case is that NewSup has only Constant terms.
	%%
	if isempty(NewSup) || ~any(any(NewSup,2))
    	%blockStruct = [];
    	%tempList = [];
	else
    	if Newm1 == 1
        	blockStruct = [blockStruct,-1];
    	else
        	blockStruct = [blockStruct,Newm1];
    	end

    	noSup = length(NewRow);
    	blockNo = repmat(nBlock+i,noSup, 1);
    	blockNo = blockNo(:);
    	idx = (1:noSup)';
    	[nonzeroIdx,temp,nonzeroVal] = find(NewCoef(:));
		idx = idx(nonzeroIdx);
    	blockNo = blockNo(nonzeroIdx);
    	tempList = [tempList;NewSup(idx,:),blockNo,NewRow(idx),NewCol(idx),nonzeroVal];
	end
end
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tempList,blockStruct] = make_varList_tC3(bSk,iEPSk,nBlock,...
    removeSup)

% This function returns information of Polynomial SDP.
% That is, monomials in this Polynomial SDP and their block number,
% row and column indices, and coefficients.
%
%
% <Return Value>
% 'tempList' has nDim+4 columns and 'blockstruct'
%
% tempList(:,1:nDim) is monomial that appear in the Polynomial SDP
% tempList(:,nDim+1) is the block No that the monomial belongs.
% tempList(:,nDim+2:nDim+3) are row and column indices in the
% matrix of block No
% tempList(:,end) is coefficients.
%
% blockStruct is the size of Polynomial SDP.
%
%
if isempty(bSk)
    tempList = [];
    blockStruct = [];
    return;
end
nDim = iEPSk.dimVar;
sCone = iEPSk.sizeCone;
m0 = size(iEPSk.supports,1);
m1 = size(bSk,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%start of finding index of upper part of moment matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Index = tril(ones(m1*sCone));
[col,row] = find(Index);
LenOneDimIdx = length(row);
sIdx = (row-1)*sCone*m1 + col;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%end of finding index of upper part of moment matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%start of Make moment matrix<--------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% finding p-th dimension support of moment matrix%%%

MomSup = sparse(LenOneDimIdx,nDim);

UsedVar = find(any(bSk,1));
idx2 = repmat((1:m1),sCone,m1*sCone);
idx1 = repmat((1:m1),m1*sCone^2,1);
idx1 = idx1(sIdx);
idx2 = idx2(:);
idx2 = idx2(sIdx);
Sup1 = bSk(idx1,UsedVar);
Sup2 = bSk(idx2,UsedVar);
MomSup(:,UsedVar) = Sup1 + Sup2;
clear Sup1 Sup2 idx1 idx2 UsedVar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%end of Make moment matrix---------->
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%start of making moment mat.\kron \F(\x)
%Output is SDPA sparse format(Only need upper part)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CoefList = [];
extIdx = repmat((1:sCone),1,m1);
for j=1:m0
    Coef = reshape(iEPSk.coef(j,:),sCone,sCone);
    Coef = Coef(extIdx,extIdx);
    Coef = Coef(:);
    Coef = Coef(sIdx);
    CoefList = [CoefList;Coef];
end
[nonZeroIdx,tmp,nonZeroVal] = find(CoefList);
idx = repmat((1:m0),LenOneDimIdx,1);
idx = idx(:);
idx = idx(nonZeroIdx);
momidx = repmat((1:LenOneDimIdx)',m0,1);
momidx = momidx(nonZeroIdx);
Sup = MomSup(momidx,:) + iEPSk.supports(idx,:);
blockNo = repmat(nBlock+1, length(nonZeroIdx), 1);
row = row(momidx);
col = col(momidx);

blockStruct = m1*sCone;
tempList = [Sup,blockNo,row,col,nonZeroVal];
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [SDPA,xIdxVec] = varListtoSDPA(varList, SDPA)

% This function converts an input 'varList', which has all
% monomials, their coefficients and indices, into SDPA sparse
% format.
%
% This function outputs 'SDPA' and 'xIdxVec'.
%
% The first is an SDPA sparse format of SDP derived from given POP.
% The second is all monomials in Polynomial SDP that sorted in
% 'grevlex' order.
%
% SDPA is the structure with the following fields:
%
% SDPA.blockStruct: the block structre vector
% SDPA.nBlock:      the number of block
% SDPA.mDim:        the nubmer of the primal variables
% SDPA.bVect:       in Primal SDP, coefficients of objective function.
% SDPA.typeCone:    indicate what the block derived from
%                   real-valued nonnegative constraint,
%                   equation or SDP constraint.
% SDPA.sparseMatrix:this is matrix with 5 columns
%                   First column is the number of matrix.
%                   Second column is the number of block.
%                   Third column is the row index.
%                   Fourth column is the column index.
%                   Last column is the coefficient.

%%
%% Assign all monomials to numbers. 'N' represents number.
%%
n = size(varList,2);
nDim = n-4;
%
% in order to sort in 'grevlex' order
%
monoList = [sum(varList(:,1:nDim),2),-varList(:,1:nDim)];
[Out, M,N] = unique(monoList,'rows');
[Row,Col,Val] = find(Out(:,2:nDim+1));
if Out(1,1) == 0
    %
    % in the case that Polynomial SDP has constant.
    %
    varNumber = (N-1)';
else
    varNumber = N';
end
xIdxVec = sparse(Row,Col,-Val,size(Out,1),nDim);
%full(xIdxVec)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Making SDPA sparse format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SDPA.mDim = max(varNumber);
bVectIndices = find(varList(:,nDim+1)==0);
bVectVarNumber = varNumber(bVectIndices);
SDPA.bVect = sparse(1,SDPA.mDim);
SDPA.bVect(1,bVectVarNumber) = varList(bVectIndices,nDim+4)';

sparseMatIndiceds = find(varList(:,nDim+1) > 0);
matrixNumber = varNumber(sparseMatIndiceds)';
blockNumber = varList(sparseMatIndiceds,nDim+1);
rowIndex = varList(sparseMatIndiceds,nDim+2);
colIndex = varList(sparseMatIndiceds,nDim+3);
value = varList(sparseMatIndiceds,nDim+4);

noOfTerms = length(matrixNumber);
SDPA.sparseMatrix = zeros(noOfTerms,5);
I = find(matrixNumber==0);
if isempty(I)
    SDPA.sparseMatrix(:,1) = matrixNumber;
    SDPA.sparseMatrix(:,2) = blockNumber;
    SDPA.sparseMatrix(:,3) = rowIndex;
    SDPA.sparseMatrix(:,4) = colIndex;
    SDPA.sparseMatrix(:,5) = value;
else
    J = find(matrixNumber);
    SDPA.sparseMatrix(1:length(I),1) = matrixNumber(I);
    SDPA.sparseMatrix(1:length(I),2) = blockNumber(I);
    SDPA.sparseMatrix(1:length(I),3) = rowIndex(I);
    SDPA.sparseMatrix(1:length(I),4) = colIndex(I);
    SDPA.sparseMatrix(1:length(I),5) = -value(I);
    SDPA.sparseMatrix(length(I)+1:end,1) = matrixNumber(J);
    SDPA.sparseMatrix(length(I)+1:end,2) = blockNumber(J);
    SDPA.sparseMatrix(length(I)+1:end,3) = rowIndex(J);
    SDPA.sparseMatrix(length(I)+1:end,4) = colIndex(J);
    SDPA.sparseMatrix(length(I)+1:end,5) = value(J);
end
clear varList;
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [NewSup,NewRow,NewCol,NewCoef] = AddSupport(MomSup,Sup,Row,Col,Coef)

% This function produces support, index and coefficient of valid
% inequality counstruced by multiplying localizing matrix and given
% constraint.
%
% MomSup, Row and Col are information for localizing matrix.
%

%%
%% Find index
%%
m0 = size(Sup,1);
noSup = size(MomSup,1);
tmpidx1 = repmat((1:noSup),1,m0);
NewRow = Row(tmpidx1');
NewCol = Col(tmpidx1');
tmpidx2 = repmat((1:m0),noSup,1);
tmpidx2 = tmpidx2(:);

%%
%% Find Coef. Coef only is extened becouse coef of localizing
%% matrix is all 1.
%%
NewCoef = Coef(tmpidx2);
%%
%% Add monomial of constraint and localizing matrix.
%%
UsedVar = find(any(Sup,1));
NewSup = MomSup(tmpidx1,:);
NewSup(:,UsedVar) = NewSup(:,UsedVar) + Sup(tmpidx2,UsedVar);


return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% M. Kojima, 02/06/2005
% inEqPolySys ---> anIneqPoly
%
function [vec,tCone] = SOCPtoSDP(anIneqPoly)

% This function converts an SOCP constraint (typeCone=2) to an SDP constraint.
%
tCone = 3;
newCoef = [];
for i=1:size(anIneqPoly.supports,1)
    sMat = anIneqPoly.coef(i,1) * speye(anIneqPoly.sizeCone);
    sMat(1,2:anIneqPoly.sizeCone) = anIneqPoly.coef(i,2:anIneqPoly.sizeCone);
    sMat(2:anIneqPoly.sizeCone,1) = anIneqPoly.coef(i,2:anIneqPoly.sizeCone)';
    newCoef = [newCoef; reshape(sMat,1,anIneqPoly.sizeCone*anIneqPoly.sizeCone)];
end
vec = newCoef;
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Sup, Row,Col] = makeMoment(basisSupports)

% This function returns monomials appeared in Moment matrix and
% their row and column indeces.
%
% Sup is monomial set.
% Row is row index.
% Col is column index.
%

%%
%% Find row and column indices.
%%
[m1,nDim] = size(basisSupports);
[Col,Row] = find(tril(ones(m1)));

Sup = basisSupports(Row,:);
%%
%% Use sparsity of basisSupports.
%%
UsedVar = find(any(basisSupports,1));
Sup(:,UsedVar) = Sup(:,UsedVar) + basisSupports(Col,UsedVar);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Sup,Row,Col,Coef,Newsize] ...
    = removeComp(Sup,Row,Col,Coef,oldsize,CompSup)
% This function removes from Polynomial SDP constrants some
% monomials that become zero by substituting complementarities.

Newsize = oldsize;
if ~isempty(CompSup)
    %%
    %% Search Removed terms
    %%
    SupSet = spones(Sup);
    SupSet = SupSet*CompSup';
    SupSet = max(SupSet,[],2);
    removeIdx = find(SupSet == 2);
    %%
    %% Remove terms
    %%
    if ~isempty(removeIdx)
        Sup(removeIdx,:) = [];
        Row(removeIdx) = [];
        Col(removeIdx) = [];
        Coef(removeIdx,:) = [];

        %%
        %% Reorder row and column indices.
        %%
        if ~isempty(Row)
            diagIdx = find(Row == Col);
            [newIdx,M,N] = unique(Row(diagIdx));
            Newsize = length(newIdx);
            if Newsize < oldsize
                [M,N,Row] = unique(Row);
                [M,N,Col] = unique(Col);
            end
        else
            Newsize = 0;
        end
    end
end


return;

function SDPA = aggregateBlocks(SDPA)

% This function gather multiple diagonal blocks derived from
% equation and inequality with size of basisSupport = 1 into one block.
%

if isfield(SDPA,'typeCone')
    eq_block_idx = find(SDPA.typeCone == -1);
    lp_block_idx = find(SDPA.typeCone ==  1);
    sdp_block_idx = find(SDPA.typeCone == 3);
else
    eq_block_idx = find(SDPA.blockStruct < -1);
    lp_block_idx = find(SDPA.blockStruct == -1);
    sdp_block_idx = find(SDPA.blockStruct > 1);
end

no_eq_block = length(eq_block_idx);
no_lp_block = length(lp_block_idx);
no_sdp_block = length(sdp_block_idx);

%%
%% renumber block No. with the order 'equation', 'linear
%% inequality', 'SDP'
%%

matNo(eq_block_idx) = (1:no_eq_block);
matNo(lp_block_idx) = no_eq_block+(1:no_lp_block);
matNo(sdp_block_idx) = no_eq_block+no_lp_block+(1: ...
    no_sdp_block);
SDPA.sparseMatrix(:,2) = matNo(SDPA.sparseMatrix(:,2));
SDPA.blockStruct = [SDPA.blockStruct(eq_block_idx),...
    SDPA.blockStruct(lp_block_idx), ...
    SDPA.blockStruct(sdp_block_idx)];

%%
%% gather some block derived from polynomial equations or
%% inequalities into one block.
%%
%% Equation Part
%%
if ~isempty(eq_block_idx)
    SDPAversion = 6.2;
    if SDPAversion <= 6.2
        %% This part gathers blocks derived from polynomial equations
        %% into one block and separate into the Positive and Negative part.
        %% The first half of one block is the Positive part.
        %% The latter half of one block is the Negative part.

        EqBlockIdx  = find(SDPA.sparseMatrix(:,2) <= no_eq_block);
        EqBlockSt_half = -SDPA.blockStruct(1:no_eq_block)'/2;
        EqBlockSt = [EqBlockSt_half;EqBlockSt_half];
        EqBlockSt = cumsum(EqBlockSt);
        EqBlockSt(end) = 0;
        EqBlockSt = circshift(EqBlockSt,1);

        EqPositiveIndex = find(SDPA.sparseMatrix(EqBlockIdx,3) <= ...
            EqBlockSt_half(SDPA.sparseMatrix(EqBlockIdx,2)));
        EqNegativeIndex = find(SDPA.sparseMatrix(EqBlockIdx,3) > ...
            EqBlockSt_half(SDPA.sparseMatrix(EqBlockIdx,2)));
        EPIdx = EqBlockIdx(EqPositiveIndex);
        ENIdx = EqBlockIdx(EqNegativeIndex);

        SDPA.sparseMatrix(EPIdx,3) =  EqBlockSt(SDPA.sparseMatrix(EPIdx,2)) ...
            +SDPA.sparseMatrix(EPIdx,3);
        SDPA.sparseMatrix(ENIdx,3) = - ...
            EqBlockSt_half(SDPA.sparseMatrix(ENIdx,2)) + ...
            EqBlockSt(SDPA.sparseMatrix(ENIdx,2)+no_eq_block) + ...
            SDPA.sparseMatrix(ENIdx,3);
        SDPA.sparseMatrix(EPIdx,4) = SDPA.sparseMatrix(EPIdx,3);
        SDPA.sparseMatrix(ENIdx,4) = SDPA.sparseMatrix(ENIdx,3);
    else
        %% In new SDPA version, we need not to separate the Positive
        %% part and Negative Part since equation part will be defined.
        %%
        EqBlockIdx  = find(SDPA.sparseMatrix(:,2) <= no_eq_block);
        EqBlockSt = -SDPA.blockStruct((1:no_eq_block))';
        EqBlockSt = cumsum(EqBlockSt);
        EqBlockSt(end) = 0;
        EqBlockSt = circshift(EqBlockSt,1);
        SDPA.sparseMatrix(EqBlockIdx,3) =  EqBlockSt(SDPA.sparseMatrix(EqBlockIdx,2)) ...
            +SDPA.sparseMatrix(EqBlockIdx,3);
        %SDPA.sparseMatrix(EqBlockIdx,4) = SDPA.sparseMatrix(EqBlockIdx,3);
    end
else
    EqBlockIdx = [];
end
%%
%% gather some block derived from polynomial equations or
%% inequalities into one block.
%%
%% Linear inequality Part
%%
if ~isempty(lp_block_idx)
    LpBlockIdx = find(SDPA.sparseMatrix(:,2) <= no_lp_block+no_eq_block);
    LpBlockIdx = setdiff(LpBlockIdx,EqBlockIdx);

    LpBlockNo = SDPA.sparseMatrix(LpBlockIdx,2)-no_eq_block;
    LpBlockSt = -SDPA.blockStruct(no_eq_block+(1:no_lp_block))';
    LpBlockSt = cumsum(LpBlockSt);
    LpBlockSt(end) = 0;
    LpBlockSt = circshift(LpBlockSt,1);
    SDPA.sparseMatrix(LpBlockIdx,3) =  LpBlockSt(LpBlockNo) ...
        +SDPA.sparseMatrix(LpBlockIdx,3);
    SDPA.sparseMatrix(LpBlockIdx,4) = ...
        SDPA.sparseMatrix(LpBlockIdx,3);
else
    LpBlockIdx = [];
end
%%
%% gather some block derived from polynomial equations or
%% inequalities into one block.
%%
%% This part has 2 procces
%% 1. Sdp Part. Here, only move block index of SDP
%% 2. Compute typeCone, nBlock and blockStruct
%%
SdpBlockIdx = find(SDPA.sparseMatrix(:,2) > no_lp_block+ ...
    no_eq_block);

if ~isempty(eq_block_idx) && ~isempty(lp_block_idx)
    SDPA.sparseMatrix(EqBlockIdx,2) = 1;
    SDPA.sparseMatrix(LpBlockIdx,2) = 2;
    SDPA.typeCone = [-1, 1, SDPA.typeCone(sdp_block_idx)];
    SDPA.sparseMatrix(SdpBlockIdx,2) = 2-no_eq_block-no_lp_block+ ...
        SDPA.sparseMatrix(SdpBlockIdx,2);
    SDPA.nBlock = 2 + no_sdp_block;
    SDPA.blockStruct = [sum(SDPA.blockStruct(1:no_eq_block)),...
        sum(SDPA.blockStruct(no_eq_block+(1:no_lp_block))), ...
        SDPA.blockStruct(no_eq_block+no_lp_block+(1: ...
        no_sdp_block))];
elseif ~isempty(eq_block_idx) && isempty(lp_block_idx)
    SDPA.sparseMatrix(EqBlockIdx,2) = 1;
    SDPA.typeCone = [-1, SDPA.typeCone(sdp_block_idx)];
    SDPA.sparseMatrix(SdpBlockIdx,2) = 1-no_eq_block+ ...
        SDPA.sparseMatrix(SdpBlockIdx,2);
    SDPA.nBlock = 1 + no_sdp_block;
    SDPA.blockStruct = [sum(SDPA.blockStruct((1:no_eq_block))), ...
        SDPA.blockStruct(no_eq_block+(1:no_sdp_block))];
elseif isempty(eq_block_idx) && ~isempty(lp_block_idx)
    SDPA.sparseMatrix(LpBlockIdx,2) = 1;
    SDPA.typeCone = [1, SDPA.typeCone(sdp_block_idx)];
    SDPA.sparseMatrix(SdpBlockIdx,2) = 1-no_lp_block+ ...
        SDPA.sparseMatrix(SdpBlockIdx,2);
    SDPA.nBlock = 1 + no_sdp_block;
    SDPA.blockStruct = [sum(SDPA.blockStruct((1:no_lp_block))), ...
        SDPA.blockStruct(no_lp_block+(1:no_sdp_block))];
end
return

% $Header: /home/waki9/CVS_DB/SparsePOPdev/subPrograms/PSDPtoLSDP.m,v 1.3 2007/01/16 14:41:14 waki9 Exp $
