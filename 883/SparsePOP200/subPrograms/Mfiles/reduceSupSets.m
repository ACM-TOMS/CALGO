%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% M. Kojima, 02/06/2006
% The main function of this module is:
%	reduceSupSets(objPoly,inEqPolySys,basisSupports);
% This module contains the following functions:
%	mulSupport(basisSup, inEqSup)
%   printGInfo(G0Key,G1Key,GcandidateKey,GstarKey,iteration)
%   printBasisInfo(noConst,basisSupports)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [basisSupports,ineqBasis] = reduceSupSets(objPoly,inEqPolySys,basisSupports,param)


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
nDim = objPoly.dimVar;
mDim = size(inEqPolySys,2);
kDim = size(basisSupports,2);
noOfinEqPolySys = size(inEqPolySys,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fe --- Even support vectors involved in the Lagrangian function

if param.boundSW == 0 && param.reduceMomentMatSW == 0
    ineqBasis = [];
    return
end
tempFe = cell(1,noOfinEqPolySys+1);
tempFe{1} = [objPoly.supports',sparse(nDim,1)];
for j=1:noOfinEqPolySys
    sup = mulSupport(basisSupports{j}, inEqPolySys{j});
    tempFe{j+1} = sup';
end
ineqBasis = [tempFe{1:noOfinEqPolySys+1}]';
ineqBasis = unique(ineqBasis, 'rows');
if param.boundSW == 1 && param.reduceMomentMatSW == 0
    return
elseif param.boundSW == 2 && param.reduceMomentMatSW == 0
	return
end
I = find(any(mod(ineqBasis,2),2) == 0);
Fe = ineqBasis(I,:);
if param.boundSW == 0 && param.reduceMomentMatSW == 1
    ineqBasis = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Attach a random key number to each support of Fe
rVector = rand(nDim,1);
FeKey = sort(Fe*rVector);
maxKey = max(FeKey);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initial G0Key{i},  GcandidateKey{i}, GstarKey{i}  (i=1,2,\ldots,mDim)
% The first element of G0{i}(j,:), Gcandidates{i}(j,:) and Gstar{i}(j,:)
% is a random keynumber
% The last element of G0{i}(j,:), Gcandidates{i}(j,:) and Gstar{i}(j,:) 
% is the original order in basisSupports{i+mDim}.
%
G0Key = cell(kDim-mDim,1);
GcandidateKey = cell(kDim-mDim,1);
GstarKey = cell(kDim-mDim,1);
G1Key = cell(kDim-mDim,1);
noOFpoints = cellfun('size',basisSupports,1);
for i=1:kDim-mDim
    G0Key{i} = [basisSupports{mDim+i} * rVector,(1:noOFpoints(mDim+i))'];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transfer \alpha from GcandidateKey{i} to GstarKey{i} such that 2*\alpha \in Fe
epsilon = 1.0e-12;
Ferow = size(FeKey,1);
G0idx = ones(Ferow,1);
noOFpoints = cellfun('size',G0Key,1);
for i=1:kDim-mDim
    value = 2*G0Key{i}(:,1);
    Feidx = ones(1,noOFpoints(i));
    FeMinusVal = FeKey(:,Feidx) -value(:,G0idx)';
    FeMinusVal = (abs(FeMinusVal) <= epsilon);
    FeMinusVal = any(FeMinusVal,1);
    OneIdx = find(FeMinusVal);
    GstarKey{i} = G0Key{i}(OneIdx,:);
    ZeroIdx = find(FeMinusVal == 0);
    I = find(value(ZeroIdx,:) - maxKey <=epsilon);
    GcandidateKey{i} = G0Key{i}(ZeroIdx(I),:);
    G1Key{i} = [GstarKey{i}; GcandidateKey{i}];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reduce elements in GcandidateKey{i} and G1Key{i} recursively until
% contSW = 0;
printSW = 0;
if printSW == 1
    iteration = 0;
    printGInfo(G0Key,G1Key,GcandidateKey,GstarKey,iteration)
end
contSW = 1;
List = cell(kDim-mDim,1);
removeIdx = cell(kDim-mDim,1);
remainIdx = cell(kDim-mDim,1);
while contSW == 1
    contSW = 0;
    % Constructing a check list
    noOFpoints = cellfun('size',G1Key,1);
    for i=1:kDim-mDim
        oneBlVal = G1Key{i}(:,ones(noOFpoints(i),1))';
        oneBlVal = oneBlVal + oneBlVal';
        [I,J,oneBlValVec] = find(tril(oneBlVal,-1));
        oneBlIdx = repmat(i,noOFpoints(i)*(noOFpoints(i)-1)/2,1);
        idx = repmat(2,noOFpoints(i),1);
        Mat = G1Key{i}(:,idx);
        [I,J,oneBlRowVec] = find(tril(Mat',-1));
        [I,J,oneBlColVec] = find(tril(Mat,-1));
        List{i} = [oneBlValVec,oneBlIdx,oneBlRowVec,oneBlColVec]';
    end
    checkList = [List{1:kDim-mDim}]';
    [V,newIdx,origIdx] = unique(checkList(:,1));
    valueList = checkList(newIdx,1);
    noOFpoints = cellfun('size',GcandidateKey,1);
    for i=1:kDim-mDim
        idx = repmat(i,noOFpoints(i),1);
        List{i} = [GcandidateKey{i},idx]';
    end
    tempList = [List{1:kDim-mDim}]';
    [V,II,JJ] = unique(tempList(:,1));
    for i=1:length(II)
        idx = find(JJ == i);
        I = find(abs(valueList -2*V(i)) <= epsilon);
        if isempty(I)
            for j=idx'
                p = tempList(j,3);
                removeIdx{p} = [removeIdx{p}, tempList(j,2)];
            end
            contSW = 1;
        else
            removeSW = 1;
            [tf,I] = ismember(I,origIdx);
            t = length(I);
            sup = sparse(t,nDim);
            p = 1;
            for k=I'
                BlIdx = checkList(k,2);
                rowIdx = checkList(k,3);
                colIdx = checkList(k,4);
                sup(p,:) = basisSupports{BlIdx+mDim}(rowIdx,:) + basisSupports{BlIdx+mDim}(colIdx,:);
                p = p+1;
            end
            
            for j = idx'
                q = tempList(j,3);
                jj = tempList(j,2);
                ridx = repmat(jj,t,1);
                sup2 = 2*basisSupports{q+mDim}(ridx,:);
                sup2 = abs(sup2 - sup);
                sup2 = (sup2 > 0);
                sup2 = any(sup2,2);
                I = find(sup2 == 0);
                if ~isempty(I)
                    removeSW = 0;
                    for k=idx'
                        p = tempList(k,3);
                        remainIdx{p} = [remainIdx{p}, tempList(k,2)];
                    end
                    break;
                end
            end
            if removeSW ~= 0
                for j=idx'
                    p = tempList(j,3);
                    removeIdx{p} = [removeIdx{p}, tempList(j,2)];
                end
                contSW = 1;
            end
        end
    end
    if contSW == 1
        contSW = 0;
        for i=1:kDim-mDim
            if ~isempty(removeIdx{i}) && isempty(remainIdx{i})
                G1Key{i} = GstarKey{i};
                GcandidateKey{i} = [];
                removeIdx{i} = [];
            elseif ~isempty(removeIdx{i}) && ~isempty(remainIdx{i})
                idx = sort(remainIdx{i});
                GcandidateKey{i} = G0Key{i}(idx,:);
                G1Key{i} = [GcandidateKey{i}; GstarKey{i}];
                contSW = 1;
                removeIdx{i} = [];
                remainIdx{i} = [];
            end
        end
    end
    if printSW == 1
        iteration = iteration + 1;
        printGInfo(G0Key,G1Key,GcandidateKey,GstarKey,iteration)
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if printSW == 1
    printBasisInfo(noOfinEqPolySys,basisSupports)
end
% <--- printing information on the basisSupports
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:kDim-mDim
    sortCol = sort(G1Key{i}(:,2));
    basisSupports{i+mDim} = basisSupports{i+mDim}(sortCol,:);
    UsedVar = find(any(basisSupports{i+mDim},1));
    [SupSet, I] = monomialSort(basisSupports{i+mDim}, UsedVar, 'grevlex');
    basisSupports{i+mDim} = basisSupports{i+mDim}(I,:);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if printSW == 1
    printBasisInfo(noOfinEqPolySys,basisSupports)
end
% <--- printing information on the basisSupports
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sup = mulSupport(basisSup, inEqPolySys)
%
% This function constructs all monomials in valid SDP constrants
% made by multiplying inequality and Moment matrix.
%
% Moreover, choose only monomials that belong in Fe (all elements
% are even.)
%

[m0,nDim0] = size(basisSup);
[m1,nDim1] = size(inEqPolySys.supports);
if nDim0 ~= nDim1
    error('dimension of sup2 is different form sup1!');
else
    nDim = nDim0;
end
usedVarL = find(any(inEqPolySys.supports,1));

if inEqPolySys.typeCone ~= -1
    %%
    %% find all monomials appeared in  Moment matrix
    %%
    mat = ones(m0);
    [col,row] = find(triu(mat));
    Msup = sparse(length(row),nDim);
    usedVarM = find(any(basisSup,1));
    Msup(:,usedVarM) = basisSup(row,usedVarM) + basisSup(col,usedVarM);
    
    %% find all monomials appeared in valid SDP
    m0 = size(Msup,1);
    idx1 = repmat((1:m0),1,m1);
    sup = Msup(idx1,:);
    idx2 = repmat((1:m1),m0,1);
    idx2 = idx2(:);
    sup(:,usedVarL) = sup(:,usedVarL) + inEqPolySys.supports(idx2,usedVarL);
elseif inEqPolySys.typeCone == -1
    %%
    %% find all monomials appeared in valid SDP
    %%
    row = repmat((1:m0),1,m1);
    sup = basisSup(row,:);
    idx = repmat((1:m1),m0,1);
    idx = idx(:);
    sup(:,usedVarL) = sup(:,usedVarL) + inEqPolySys.supports(idx,usedVarL);
end
return

function printBasisInfo(noConst,basisSupports)
fprintf('** basisSupports\n');
cc = size(basisSupports,2);
for p=noConst+1:cc
    rowSize = size(basisSupports{p},1);
    fprintf('%2d : \n',p);
    for i=1:rowSize
        fprintf('    ');
        vec = basisSupports{p}(i,:);
        if issparse(vec)
            fprintf('%2d',full(vec'));
        else
            fprintf('%2d',vec');
        end
        fprintf('\n');
    end
end
return

function printGInfo(G0Key,G1Key,GcandidateKey,GstarKey,iteration)
fprintf('* %3d *\n',iteration);
for i=1:size(G0Key,1)
    fprintf('%2d: ',i);
    fprintf('#G0 = %3d, ',size(G0Key{i},1));
    fprintf('#G1 = %3d, ',size(G1Key{i},1));
    fprintf('#Gcandidate = %3d, ',size(GcandidateKey{i},1));
    fprintf('#Gstar = %3d\n',size(GstarKey{i},1));
end
fprintf('\n');
return

% $Header: /home/waki9/CVS_DB/SparsePOPdev/subPrograms/reduceSupSets.m,v 1.2 2007/01/16 08:17:56 waki9 Exp $
