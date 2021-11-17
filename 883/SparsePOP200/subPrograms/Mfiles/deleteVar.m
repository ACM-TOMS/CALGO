function [NewobjPoly,NewineqPolySys,Newlbd,Newubd,fixedVar] = deleteVar(objPoly,ineqPolySys,lbd,ubd,param)

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
%t = cputime;
%{
disp('objPoly'); 
full([objPoly.supports,objPoly.coef])
disp('ineqPolySys');
for i=1:size(ineqPolySys,2)
	full([ineqPolySys{i}.supports,ineqPolySys{i}.coef])
end
[lbd;ubd]
%}


count = 1;
nDim = objPoly.dimVar;
fixedVar = [];
aliveIdx = 1:nDim;
while count == 1
    [NewobjPoly,NewineqPolySys,Newlbd,Newubd,NewaliveIdx,NewfixedVar] = deleteVarMain(aliveIdx,objPoly,ineqPolySys,lbd,ubd,param);
    %%
	%%
	%% display for debug
	%%
	%{
	disp('objPoly'); 
	full([NewobjPoly.supports,NewobjPoly.coef])
    disp('ineqPolySys');
	for i=1:size(NewineqPolySys,2)
    	full([NewineqPolySys{i}.supports,NewineqPolySys{i}.coef])
    end
	[Newlbd;Newubd]
	%}
	%%
    if isempty(NewfixedVar)
        count = 0;
    else
	aliveIdx = NewaliveIdx;
        fixedVar = [fixedVar;NewfixedVar];
        objPoly = NewobjPoly;
        ineqPolySys = NewineqPolySys;
        lbd = Newlbd;
        ubd = Newubd;
    end
end
if ~isempty(fixedVar)
    aliveIdx = ones(1,nDim);
    aliveIdx(fixedVar(:,1)) = 0;
% Kojima 06/12/2007 ---> 
%	showDeletedVar(aliveIdx);
% 	if NewobjPoly.dimVar == 0
% 		fprintf('\n');
% 		fprintf('## All variables x_i in POP are removed.            ##\n');
% 		fprintf('## SparsePOP does not apply SDP relaxtion into POP. ##\n');
% 	end
% Kojima <--- 06/12/2007
end
if isempty(fixedVar)
	NewobjPoly = objPoly;
	t = size(ineqPolySys,2);
	NewineqPolySys = cell(1,t);
	for i=1:size(ineqPolySys,2)
		NewineqPolySys{i} = ineqPolySys{i};
	end
	Newlbd = lbd;
	Newubd = ubd;
end



%t = cputime - t;
%fprintf('Cpu Time of deleteVar = %3.3e\n',t);

return


function [NewobjPoly,NewineqPolySys,Newlbd,Newubd,NewaliveIdx,fixedVar] = deleteVarMain(aliveIdx,objPoly,ineqPolySys,lbd,ubd,param)


nDim = objPoly.dimVar;
fixedVar = [];
lbdIdx = find(lbd > -1.0e+8);
ubdIdx = find(ubd <  1.0e+8);
Idx = intersect(lbdIdx,ubdIdx);

%%% Step 1: Get variables whose lower bound is same as the upper bounds.
%%%         And, substituting into objective and constraint functions

if ~isempty(Idx)
    d = ubd(Idx)-lbd(Idx);
    zeroIdx = find(abs(d) <  1.0e-8);
	negaIdx = find(d < -1.0e-8);
	if ~isempty(negaIdx)
        error('ubd < lbd at some index.');
    elseif ~isempty(zeroIdx)
        Idx = Idx(zeroIdx);
        val = lbd(zeroIdx);
        [aliveIdx,fixedVar,varval] = getDeletedVar(Idx,val,aliveIdx,fixedVar);
        [NewobjPoly,NewineqPolySys,Newlbd,Newubd] = SubstituteIntoProblems(objPoly,ineqPolySys,lbd,ubd,fixedVar); 
        %dips('step1');
		%showDeletedVar(aliveIdx);
    else
        NewobjPoly = objPoly;
        NewineqPolySys = ineqPolySys;
        Newubd = ubd;
        Newlbd = lbd;
    end
else
    NewobjPoly = objPoly;
    NewineqPolySys = ineqPolySys;
    Newubd = ubd;
    Newlbd = lbd;
end
%disp('step1');
%fixedVar
%%% Step 2: Eliminate redundant complementarity constraints %%%
if 1 == 1
    Idx = getComplimentarity(NewineqPolySys,Newlbd,Newubd,param);
	if ~isempty(Idx)
        val = zeros(1, length(Idx));
        [aliveIdx,fixedVar,varval] = getDeletedVar(Idx,val,aliveIdx,fixedVar);
		[NewobjPoly,NewineqPolySys,Newlbd,Newubd] = SubstituteIntoProblems(NewobjPoly,NewineqPolySys,Newlbd,Newubd,varval);
        %disp('step2');
		%showDeletedVar(aliveIdx)
    end
end
%disp('step2');
%fixedVar
%%% Step 3: Getting "x_i = const", and substituting
if 1 == 1
    zeroIdx = [];
    NonzeroIdx = [];
    for i=1:size(NewineqPolySys,2)
        typeCone = NewineqPolySys{i}.typeCone;
        sizeCone = NewineqPolySys{i}.sizeCone;
		degree   = NewineqPolySys{i}.degree;
        if typeCone == -1 && sizeCone == 1 && degree == 1 
            noTerms = NewineqPolySys{i}.noTerms;
            if noTerms == 1 %% x_i = 0
                varidx = find(NewineqPolySys{i}.supports);
                if length(varidx) == 1
                    zeroIdx = [zeroIdx;varidx,0];
                end
            elseif noTerms== 2 %% x_i = nonzero
                [row,col] = find(NewineqPolySys{i}.supports);
                FstSup = NewineqPolySys{i}.supports(1,:);
                SndSup = NewineqPolySys{i}.supports(2,:);
                if ~any(FstSup,2) && length(row) == 1 
                    %% FstSup = zeros
                    varidx = find(SndSup);
                    coef = -NewineqPolySys{i}.coef(1,1)/NewineqPolySys{i}.coef(2,1);
                    NonzeroIdx = [NonzeroIdx;varidx,coef];
                elseif ~any(SndSup) && length(row) == 1
                    varidx = find(FstSup);
                    coef = -NewineqPolySys{i}.coef(2,1)/NewineqPolySys{i}.coef(1,1);
                    NonzeroIdx = [NonzeroIdx;varidx,coef];                     
                end
            end
        end
    end
    Idxval = [zeroIdx;NonzeroIdx];
    if ~isempty(Idxval)
        idx = Idxval(:,1);
        val = Idxval(:,2);
        [aliveIdx,fixedVar,varval] = getDeletedVar(idx,val,aliveIdx,fixedVar);
		[NewobjPoly,NewineqPolySys,Newlbd,Newubd] = SubstituteIntoProblems(NewobjPoly,NewineqPolySys,Newlbd,Newubd,varval);
	%disp('step3');
	%showDeletedVar(aliveIdx);
    end
end
NewaliveIdx = aliveIdx;
%disp('step3');
%fixedVar
return


function [NewobjPoly,NewineqPolySys,Newlbd,Newubd] = SubstituteIntoProblems(objPoly,ineqPolySys,lbd,ubd,fixedVar)

zeroIdx = fixedVar(:,1);
nDim = objPoly.dimVar;
remainIdx = setdiff([1:nDim], zeroIdx);
val = fixedVar(:,2);
count = 0;

%% objPoly
NewobjPoly = objPoly;
mono.typeCone = 1;
mono.sizeCone = 1;
mono.dimVar = length(zeroIdx);
mono.noTerms = 1;
for i=1:objPoly.noTerms
    mono.supports = objPoly.supports(i,zeroIdx);
    mono.coef = objPoly.coef(i,:);
    if size(mono.supports,1) == 1&& any(mono.supports,2) == 0
        value = mono.coef;
    else
        [value,temp] = evalPolynomials(mono,val);
    end
    NewobjPoly.coef(i,:) = value;
end
NewobjPoly.dimVar = objPoly.dimVar - length(zeroIdx);
NewobjPoly.supports = objPoly.supports(:,remainIdx);
NewobjPoly = simplifyPolynomial(NewobjPoly);

%% ineqPolySys
for i=1:size(ineqPolySys,2)
    mono.sizeCone = ineqPolySys{i}.sizeCone;
    for j=1:ineqPolySys{i}.noTerms
        mono.supports = ineqPolySys{i}.supports(j,zeroIdx);
        mono.coef = ineqPolySys{i}.coef(j,:);
        if size(mono.supports,1) == 1&& any(mono.supports,2) == 0
            value = mono.coef;
        else
            [value,temp] = evalPolynomials(mono,val);
        end
        ineqPolySys{i}.coef(j,:) = value;
    end
    ineqPolySys{i}.dimVar = ineqPolySys{i}.dimVar - length(zeroIdx);
    ineqPolySys{i}.supports = ineqPolySys{i}.supports(:,remainIdx);
    ineqPolySys{i} = simplifyPolynomial(ineqPolySys{i});

	if ~isempty(ineqPolySys{i}) 
		if  ~isempty(ineqPolySys{i}.supports)
	    	if ineqPolySys{i}.noTerms  > 1
    	   		count = count + 1;
	    	    NewineqPolySys{count} = ineqPolySys{i};
	    	elseif ineqPolySys{i}.noTerms == 1
    	    	if any(ineqPolySys{i}.supports,2) == 0
					I = [];
					for k=1:ineqPolySys{i}.sizeCone
 						if abs(ineqPolySys{i}.coef(1,k)) > 1.0e-8
							I = [I, k];		
                        end
                    end
                    if ~isempty(I)
        	    		count = count + 1;
            			NewineqPolySys{count} = ineqPolySys{i};
						NewineqPolySys{count}.sizeCone = length(I);
						NewineqPolySys{count}.coef = ineqPolySys{i}.coef(:,I);
					end
                else
                    count = count + 1;
                    NewineqPolySys{count} = ineqPolySys{i};
                end
			end
		end
    end

end
if count == 0
	NewineqPolySys = cell(1,0);
end
%for i=1:size(NewineqPolySys,2)
%       full([NewineqPolySys{i}.supports, NewineqPolySys{i}.coef]) 
%end


%% ubd
Newubd = ubd(remainIdx);
%% lbd
Newlbd = lbd(remainIdx);

return

function Idx = getComplimentarity(NewineqPolySys,Newlbd,Newubd,param)

Idx = [];
%%% Step 1 --- Pick up complementarity constraints %%%
param.complementaritySW = 1;
[CompSup,ConstraintInfo] = separateSpecialMonomial(NewineqPolySys,param);
if isempty(CompSup)
    return;
end
%%% Step 2 --- Pick up supports of all comprimentarity constraints %%%

%%% Step 3 --- Pick up variables whose box constraint does not have 0 %%%

vec = (Newlbd.*Newubd > 1.0e-8);
if ~any(vec,2)
    return;
end

%%% Step 4 --- Find complimentarity constraints with positive or negative
%%% variables %%%
one = find(CompSup*vec' == 1);
if ~ isempty(one)
    [row,col] = find(CompSup(one,:)-repmat(vec,length(one),1) > 0);
    Idx = unique(col);
end

two = find(CompSup*vec' == 2);
if ~isempty(two)
    error('This problem is infeasible. Because this has complimentarity constraint %d, but variables used in this constraint can not become zero',two(1));
end
return

function [aliveIdx,fixedVar,varval] = getDeletedVar(idx,val,aliveIdx,fixedVar)
aIdx = find(aliveIdx);

if size(aIdx(idx),1) ~= 1
    aIdx(idx) = aIdx(idx)';
end
if size(val,1) ~= 1
    val = val';
end

fixedVar = [fixedVar; aIdx(idx)',val'];
aliveIdx(aIdx(idx)) = 0;

if size(idx,1) >1
   idx = idx'; 
end
varval = [idx',val'];
return

function showDeletedVar(aliveIdx)
Idx = find(aliveIdx == 0);
if ~isempty(Idx)
    fprintf('**********************************');
    fprintf('*******************************\n');
    fprintf('SparsePOP automatically finds the constraint "x(j) = constant"\n');
    fprintf('and removes the following variables from the original POP\n');
	fprintf('by substituting them into the POP.\n');

    for i=Idx
        if mod(i,10) == 1 && i ~= 11
            fprintf('%2dst variable Eliminated\n',i);
        elseif mod(i,10) == 2 && i~= 12
            fprintf('%2dnd variable Eliminated\n',i);
        elseif mod(i,10) == 3 && i~= 13
            fprintf('%2drd variable Eliminated\n',i);
        else
            fprintf('%2dth variable Eliminated\n',i);
        end
    end
    fprintf('**********************************');
    fprintf('*******************************\n');
end
return

% $Header:
