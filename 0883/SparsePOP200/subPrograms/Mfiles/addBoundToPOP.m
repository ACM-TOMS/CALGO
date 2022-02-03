function [NewineqPolySys,NewbasisSupports,boundList] = addBoundToPOP(ineqPolySys,...
    basisSupports,lbd,ubd,momentMono,ineqMono,lbdIdx,ubdIdx,param)

% This function adds lower and upper bounds for variables which appear in
% moment problem.
%
% If lower and upper bounds of all variables in POP are 0 and 1, this
% function ommits redundant bounds of variables in moment problem.
%

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
if param.boundSW == 0
   NewineqPolySys = ineqPolySys;
   NewbasisSupports = basisSupports;
   boundList = [];
   return;
end
nDim = length(ubd);
if nDim ~= length(lbd)
    error('the size of ubd is different from the size of lbd.');
end
pointer = size(ineqPolySys,2);
noLocal = pointer;
noMoment= size(basisSupports,2) - noLocal;


I = find((lbd < 1.0e-8) & (lbd > -1.0e-8));
J = find((ubd <  1.0+1.0e-8) & (ubd >  1.0-1.0e-8));

removeSW = 0;
if length(I) == nDim && length(J) == nDim 
    removeSW = 1;
end

if removeSW == 1
    [OneSup,ZeroSup] = addZeroOnebox(ineqPolySys,basisSupports, momentMono,ineqMono,lbdIdx,ubdIdx,param);
	[OneSup,J] = monomialSort(OneSup,1:nDim,'grevlex');
	[ZeroSup,J] = monomialSort(ZeroSup,1:nDim,'grevlex');
    t1 = size(OneSup,1);
    t2 = size(ZeroSup,1);
	%full(OneSup)
	%full(ZeroSup)
	if t1 == 0 && t2 == 0
        boundList = [];
    else
		temp = sparse(2*t1,nDim);
		I = 2*(1:t1);
		temp(I,:) = OneSup;
        boundList = sparse(2*t1+t2,nDim+4);
        boundList(:,1:nDim) =[ZeroSup;temp];
        ridx = t2+(1:t1);
		ridx1 = repmat(ridx,2,1);
		ridx1 = ridx1(:);
        ridx2 = (1:t2)';
        boundList(:,nDim+2) = [ridx2;ridx1];
        boundList(:,nDim+3) = [ridx2;ridx1];
		coef = [1;-1];
		coef = repmat(coef,t1,1);
        boundList(:,end) = [ones(t2,1);coef];
		%{
        boundList = sparse(2*t1+t2,nDim+4);
        boundList(:,1:nDim) =[ZeroSup;sparse(t1,nDim);OneSup];
        ridx = t2+(1:t1)';
        ridx1 = [ridx;ridx];
        ridx2 = (1:t2)';
        boundList(:,nDim+2) = [ridx2;ridx1];
        boundList(:,nDim+3) = [ridx2;ridx1];
        boundList(:,end) = [ones(t2,1);ones(t1,1);-ones(t1,1)];
		%}
		%{
        boundList = sparse(2*t1+t2,nDim+4);
        boundList(:,1:nDim) =[OneSup;sparse(t1,nDim);ZeroSup];
        ridx = (1:t1)';
        ridx1 = [ridx;ridx];
        ridx2 = t1+(1:t2)';
        boundList(:,nDim+2) = [ridx1;ridx2];
        boundList(:,nDim+3) = [ridx1;ridx2];
        boundList(:,end) = [-ones(t1,1);ones(t1+t2,1)];
    	%}
	end
    NewbasisSupports = basisSupports;
    NewineqPolySys = ineqPolySys;
elseif removeSW == 0
    boundList = [];
    allMono = [momentMono;ineqMono];
    allMono = unique(allMono,'rows');
    NonLinIdx = find(sum(allMono,2) > 1);
    allMono = allMono(NonLinIdx,:);
	[allMono,J] = monomialSort(allMono,1:nDim,'grevlex');
    noMono = size(allMono,1);
	%full(allMono)
	%temp = cell(1,noMono);
    temp = cell(1,0);
	for i=1:noMono
        mono = allMono(i,:);
        I = find(mono);
        temp{i}.typeCone = 1;
        temp{i}.dimVar   = nDim;
        temp{i}.noTerms  = 2;
        temp{i}.supports = [sparse(1,nDim); mono];
        temp{i}.degree   = sum(mono(I),2);

        neglbdidx = find(lbd(I) < -1.0e-8);
        finiteubdidx = find(ubd(I)>= 1.0e+8);
        zerlbdidx = find(abs(lbd(I)) < 1.0e-8);
        % attach a lower bound to a y variable
        if isempty(neglbdidx) && isempty(finiteubdidx)
            %% To avoid n-th power of 0.
            if ~isempty(zerlbdidx)
                yLowerBound = 0;
            else
                yLowerBound = prod(power(lbd(I),mono(I)));
            end
            yUpperBound = prod(power(ubd(I),mono(I)));
            temp{i}.sizeCone = 2;
            temp{i}.coef     = [-yLowerBound yUpperBound; 1 -1];
            pointer = pointer + 1;
        elseif isempty(neglbdidx)
            %% To avoid n-th power of 0.
            if ~isempty(zerlbdidx)
                yLowerBound = 0;
            else
                yLowerBound = prod(power(lbd(I),mono(I)));
            end
            temp{i}.sizeCone = 1;
            temp{i}.coef     = [-yLowerBound; 1];
            pointer = pointer + 1;
        else
            absBd = max(abs(lbd(I)),abs(ubd(I)));
            InfubdIdx = find(absBd > 1.0e+8);
            if isempty(InfubdIdx) && any(mod(mono,2),2) == 0
                yBound = prod(power(absBd,mono(I)));
                temp{i}.sizeCone = 2;
                temp{i}.coef     = [yBound 0;-1 1];
                pointer = pointer + 1;
            elseif isempty(InfubdIdx) && any(mod(mono,2),2) == 1
                yBound = prod(power(absBd,mono(I)));
                temp{i}.sizeCone = 2;
                temp{i}.coef     = [yBound yBound;-1 1];
                pointer = pointer + 1;
            	%%full([temp{i}.supports,temp{i}.coef])
			else
                temp{i}.coef = [];
            end
        end
    end
    NewbasisSupports = cell(1,pointer+noMoment);
    NewineqPolySys = cell(1,pointer);
    for i=1:noLocal
        NewbasisSupports{i}  = basisSupports{i};
        NewineqPolySys{i} = ineqPolySys{i};
    end
    if ~isempty(temp)
        q = noLocal;
        for i=1:size(temp,2)
            if ~isempty(temp{i}) && ~isempty(temp{i}.coef) 
                if size(temp{i}.supports,1) ~= 1
					q = q+1;
                	NewineqPolySys{q} = temp{i};
                	NewbasisSupports{q}  = sparse(1,nDim);
                    %full([temp{i}.supports,temp{i}.coef])
				elseif size(temp{i}.supports,1) == 1 && any(temp{i}.supports) == 1
					q = q+1;
                	NewineqPolySys{q} = temp{i};
                	NewbasisSupports{q}  = sparse(1,nDim);
                    %full([temp{i}.supports,temp{i}.coef])
				end
            end
        end
    end
    for i=pointer+1:pointer+noMoment
        NewbasisSupports{i} = basisSupports{i-pointer+noLocal};
    end

end
return

function [OneMono,ZeroMono] = addZeroOnebox(ineqPolySys,basisSupports,momentMono,ineqMono,lbdIdx,ubdIdx,param)

% This function assumes that all variables in POP lie in [0,1].
%
% outputs:
% OneMono has variables whose upper bounds are 1.
% ZeroMono has variables whose lower bounds are 0

nDim = size(momentMono,2);
noineq = size(ineqPolySys,2);
noMoment = size(basisSupports,2);
allMono = [momentMono;ineqMono];
I = find(sum(allMono,2) <= 1);
allMono(I,:) = [];
if param.boundSW == 1
    allMono = unique(allMono,'rows');
    ZeroMono = allMono;
    OneMono = allMono;
    return
end
DeletedSup = [];
%%
%% for ZeroMono
%%
for i= 1:nDim
    if ~isempty(lbdIdx{i})
        for j=1:length(lbdIdx{i})
            sup = 2*basisSupports{lbdIdx{i}(j)};
            sup(:,i) = sup(:,i) + 1;
            DeletedSup = [DeletedSup;sup];
        end
    end
end
for i=noineq+1:noMoment
    DeletedSup =[DeletedSup;2*basisSupports{i}];
end

if isempty(DeletedSup)
    ZeroMono = allMono;
elseif ~isempty(DeletedSup) && isempty(allMono)
    ZeroMono = [];
else
    ZeroMono = setdiff(allMono,DeletedSup,'rows');
end


%%
%% for OneMono
%%
DeletedSup = [];
for i=noineq+1:noMoment
    nob = size(basisSupports{i},1);
    [I,J] = find(triu(ones(nob),1));
    basis = basisSupports{i}(I,:) + basisSupports{i}(J,:);
    DeletedSup = [DeletedSup;basis];
end
for i= 1:nDim
    if ~isempty(ubdIdx{i})
        for j=1:length(ubdIdx{i})
            sup = 2*basisSupports{ubdIdx{i}(j)};
            sup(:,i) = sup(:,i) + 1;
            DeletedSup = [DeletedSup;sup];
        end
    end
end

if ~isempty(DeletedSup) && ~isempty(allMono)
    OneMono = setdiff(allMono,DeletedSup,'rows');
elseif isempty(DeletedSup)
    OneMono = allMono;
else
    OneMono = [];
end

    
%print(OneMono,ZeroMono);

return

function print(OneMono, ZeroMono)
fprintf('# of bounds that added into Poly.SDP\n');
fprintf('x^a le 1 --- %d\n',size(OneMono,1));
if size(OneMono,2) < 100
    allWrite = 1;
else
    allWrite = 0;
end
if allWrite == 1
    for i=1:size(OneMono,1)
        for j=1:size(OneMono,2)
            fprintf('%d ',OneMono(i,j));
        end
        fprintf('\n');
    end
    fprintf('\n');
end
fprintf('x^a ge 0 --- %d\n',size(ZeroMono,1));
if allWrite == 1
    for i=1:size(ZeroMono,1)
        for j=1:size(ZeroMono,2)
            fprintf('%d ',ZeroMono(i,j));
        end
        fprintf('\n');
    end
end
return



% $Header:
