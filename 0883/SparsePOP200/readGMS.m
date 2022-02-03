function [objPoly,ineqPolySys,lbd,ubd] = readGMS(fileName,symbolicMath)
%
% readGMS
% converts GAMS scalar format into SparsePOP format
%
% Usage:
% [objPoly,ineqPolySys,lbd,ubd] = readGMS(fileName,symbolicMath);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fileName: the file name written in GAMS scalar format.
% symbolicMath: 1 if you have the Symbolic Math Toolbox provieded by
%      MathWorks. With this option, parentheses in the file are
%      expanded automatically. Default value is 0.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% objPoly, inEqPolySys, lbd, ubd form the SparsePOP format.
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
if nargin < 2
    symbolicMath = 0;
end

eqTo2ineqSW = 0;
eqTolerance = 0.0;

% Input %%%%%
% filename: GAMS data file name with .gms, for example "ex2_1_2.gms".
% eqTo2ineqSW
%   = 0;    to keep equalities as they are.
%   = 1;    to convert an equality f(x) = 0 into f(x) >= 0
%           and f(x) <=  eqTolerance.
% eqTolerance = 1e-4;
%%%%%%%%%%%%%
% Restriction:
%   1.  A line starting with '*' in the first column is regarded as a comment
%       line.
%   2.  At most one item of "Variables", "Positive Variables", "Equations",
%       constraints and bounds is contained in one line.
%	The end of one item is marked with ';'.
%	One item can be written in more than one line;
%       The first letter of each line can not be '*'.
%       For example,
%           Positive Variables x1,x2,x3,x4,
%               x5,x6;
%           e1..  +0.5*x1*x1 +0.5*x2*x2 +0.5*x3*x3 +0.5*x4*x4
%               +0.5*x5*x5 +10.5*x1 +7.5*x2
%               +3.5*x3 +2.5*x4 +1.5*x5 +10*x6 + objvar =E= 0;
%           x1.up = 1;
%       are allowed. But
%           x1.lo = 1; x1.up = 2;
%           e1..  +0.5*x1*x1 +0.5*x2*x2 +0.5*x3*x3 +0.5*x4*x4 +3.5
%                 *x3 +2.5*x4 +1.5*x5 +10*x6 + objvar =E= 0;
%       are not allowed.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output %%%%%
% objPoly
% ineqPolySys
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Variables
% noOfVariables <--- "Variables" in the GAMS file;
% varNames{k} (k=1,2,?ldots,noOfVariables)
%   <--- "Variables" in the GAMS file;
% posVarNames (k=1,2,?ldots,noOfVariables)
%   <--- "Positive Variables" in the GAMS file;
% noOfEquations
%   <--- "Equations" in the GAMS file;
% equationNames (i=1,2,?ldots,noOfEquations);
%   <--- "Equations" in the GAMS file;
% lbd(1,k) (k=1,2,?ldots,noOfVariables);
% ubd(1,k) (k=1,2,?ldots,noOfVariables);
% listOfTerms{i} (i=1,2,?ldots,noOfEquations);
%   --- the list of monomials of each equation.

ineqPolySys = [];

% Reading the GAMs file "fileName"
fileIDX = fopen(fileName,'r');
statusSW = 1;
noOfEquations = 0;
pp = 1;
posVarNames = [];
while statusSW == 1
    [statusSW,oneLine] = getOneLine(fileIDX);
    if statusSW == 1
        %		fprintf('%s?n',oneLine);
        [keyword,oneLine] = strtok(oneLine);
        %		fprintf('%s = %s?n',keyword,remLine);
        if strcmp('Variables',keyword)
            %			fprintf('?n%s = %s?n?n',keyword,remLine);
            varNames = [];
            p = 0;
            [varNames,p,moreSW] = getListOfNames(oneLine,varNames,p);
            while moreSW == 1
                [statusSW,oneLine] = getOneLine(fileIDX);
                [varNames,p,moreSW] = getListOfNames(oneLine,varNames,p);
            end
            %			varNames
            noOfVariables = size(varNames,2);
            lbd = -1.0e10* ones(1,noOfVariables);
            ubd = 1.0e10* ones(1,noOfVariables);
            fixed = lbd;
        elseif strcmp('Positive',keyword)
            [keyword,oneLine] = strtok(oneLine);
            if strcmp('Variables',keyword)
                %				posVarNames = [];
                p = 0;
                [posVarNames,p,moreSW] = getListOfNames(oneLine,posVarNames,p);
                while moreSW == 1
                    [statusSW,oneLine] = getOneLine(fileIDX);
                    [posVarNames,p,moreSW] = getListOfNames(oneLine,posVarNames,p);
                end
                %				posVarNames
            end
        elseif strcmp('Equations',keyword)
            equationNames  = [];
            p = 0;
            [equationNames,p,moreSW] = getListOfNames(oneLine,equationNames,p);
            while moreSW == 1
                [statusSW,oneLine] = getOneLine(fileIDX);
                [equationNames,p,moreSW] = getListOfNames(oneLine,equationNames,p);
            end
            %			equationNames
            noOfEquations = size(equationNames,2);
            listOfEquations = [];
        elseif pp <= noOfEquations
            if strcmp(strcat(equationNames{pp},'..'),keyword)
                if symbolicMath == 1
                    if ~isempty(strfind(oneLine,'('))
                        x1 = sym('x1');
                        loca = strfind(oneLine,'objvar');
                        if isempty(loca)
                            loca = strfind(oneLine,'=');
                            loca = loca -1;
                        else
                            loca = loca -3; % for objvar
                        end
                        oneLinetmp = char(vpa(expand(collect(oneLine(1:loca),x1)),20));
                        oneLine = strcat(oneLinetmp,oneLine(loca+1:size(oneLine,2)));
		    end
                else
			if ~isempty(strfind(oneLine,'('))
				error('Please expand parenthesises by your hand.');
			end
		end
                oneLinetmp = oneLine;		% to remove blank around *
                while ~isempty(strfind(oneLinetmp,' *')) || ~isempty(strfind(oneLinetmp,'* '))
                    if ~isempty(strfind(oneLinetmp, ' *'))
                        loca = strfind(oneLinetmp,' *');
                        loca = loca -1;
                        oneLinetmp=strcat(oneLinetmp(1:loca),oneLinetmp(loca+2:size(oneLinetmp,2)));
                    elseif ~isempty(strfind(oneLinetmp, '* '))
                        loca = strfind(oneLinetmp,'* ');
                        oneLinetmp=strcat(oneLinetmp(1:loca),oneLinetmp(loca+2:size(oneLinetmp,2)));
                    end
                end
                oneLine = oneLinetmp;
                %				oneLine
		listOfEquations{pp} = oneLine;
                pp = pp+1;
            end
        elseif (0 < noOfEquations) && (noOfEquations < pp)
            [oneVarName,bound] = strtok(keyword,'. ');
            %			oneVarName
            %			bound
            for i=1:noOfVariables
                if strcmp(oneVarName,varNames{i})
                    % 					oneLine
                    asciiVal = strtok(oneLine,' =;');
                    % 					asciiVal
                    if strcmp(bound,'.lo')
                        lbd(1,i) = str2double(asciiVal);
                    elseif strcmp(bound,'.up')
                        ubd(1,i) = str2double(asciiVal);
                    elseif strcmp(bound,'.fx')
                        fixed(1,i) = str2double(asciiVal);
                        lbd(1,i) = str2double(asciiVal);
                        ubd(1,i) = str2double(asciiVal);
                    end
                end
            end
        end
    end
end

printSW = 0;
if printSW == 1
    nnn = size(lbd,2);
    for i=1:nnn
        fprintf('%+6.2e ',lbd(1,i));
    end
    fprintf('\n');
    for i=1:nnn
        fprintf('%+6.2e ',ubd(1,i));
    end
    fprintf('\n');
end

noOfVariables = size(varNames,2);

% finding the objective row
i = 1;
temp = [];
while (i <= noOfEquations) && (isempty(temp))
    temp = findstr(listOfEquations{i},'objvar');
    if isempty(temp) ~= 1
        objRow = i;
    end
    i = i+1;
end

% analyzing list of constraints.
% separating each line of polynomial equation into the list of monomials.
listOfTerms = cell(1,noOfEquations);
eqOrIneq = cell(1,noOfEquations);
rightValue = cell(1,noOfEquations);
for i=1:noOfEquations
    [listOfTerms{i},eqOrIneq{i},rightValue{i}] = separate(listOfEquations{i});
    %% check rightValue{i}
    %% the rightValue{i} must be a constant value.
    if ~isnumeric(rightValue{i}) || isempty(rightValue{i}) || isnan(rightValue{i})
        error('The right-hand value of objective function and constraints must be a constant value.')
    end
end


% finding the objective term in the objective row
noOfTerms = size(listOfTerms{objRow},2);
p = 1;
temp = [];
while (p <= noOfTerms) && isempty(temp)
    temp = findstr(listOfTerms{objRow}{p},'objvar');
    if isempty(temp) ~= 1
        objTerm = p;
    end
    p = p+1;
end

% eliminating objvar from varNames
p=0;
for i=1:noOfVariables
    if strcmp('objvar',varNames{i}) ~= 1
        p = p+1;
        varNames{p} = varNames{i};
        lbd(p) = lbd(i);
        ubd(p) = ubd(i);
        fixed(p) = fixed(i);
    end
end
noOfVariables = noOfVariables -1;
lbd = lbd(1:noOfVariables);
ubd = ubd(1:noOfVariables);
fixed = fixed(1:noOfVariables);

printSW = 0;
if printSW == 1
    nnn = size(lbd,2);
    for i=1:nnn
        fprintf('%+6.2e ',lbd(1,i));
    end
    fprintf('\n');
    for i=1:nnn
        fprintf('%+6.2e ',ubd(1,i));
    end
    fprintf('\n');
end

q = 0;
if listOfTerms{objRow}{objTerm}(1) == '-'
    objConstant = -rightValue{objRow};
    objFunction = cell(1,noOfTerms);
    for p=1:noOfTerms
        if p ~= objTerm
            q = q+1;
            objFunction{q} = listOfTerms{objRow}{p};
        end
    end
else
    objConstant = rightValue{objRow};
    for p=1:noOfTerms
        if p ~= objTerm
            q = q+1;
            ll = length(listOfTerms{objRow}{p});
            if listOfTerms{objRow}{p}(1) == '-'
                objFunction{q} = strcat('+',listOfTerms{objRow}{p}(2:ll));
            else
                objFunction{q} = strcat('-',listOfTerms{objRow}{p}(2:ll));
            end
        end
    end
end

% eliminating the objective row from listOfTerms
q = 0;
for p=1:noOfEquations
    if p ~=objRow
        q = q+1;
        listOfTerms{q} = listOfTerms{p};
        eqOrIneq{q} = eqOrIneq{p};
        rightValue{q} = rightValue{p};
    end
end
noOfEquations = noOfEquations - 1;

debug = 0;
if debug
    fprintf('varNames:   ')
    for i=1:noOfVariables
        fprintf('%5s     ',varNames{i});
    end
    fprintf('\n');
    if isempty(posVarNames) ~= 1
        ll = size(posVarNames,2);
        fprintf('posVarNames:')
        for i=1:ll
            fprintf('%5s     ',posVarNames{i});
        end
        fprintf('\n');
    end
    fprintf('lbd   : ');
    for i=1:noOfVariables
        fprintf('%+7.2e ',lbd(i));
    end
    fprintf('\n');
    fprintf('ubd   : ');
    for i=1:noOfVariables
        fprintf('%+7.2e ',ubd(i));
    end
    fprintf('\n');
    fprintf('objFunction : ');
    ll = size(objFunction,2);
    for j=1:ll
        fprintf('%s ',objFunction{j});
    end
    fprintf('\n');
    for i=1:noOfEquations
        fprintf('%2d : ',i);
        ll = size(listOfTerms{i},2);
        for j=1:ll
            fprintf('%s ',listOfTerms{i}{j});
        end
        fprintf(' %s ',eqOrIneq{i});
        fprintf(' %+8.3e \n',rightValue{i});
    end
end

objPoly1 = convToPolynomial(noOfVariables,varNames,objFunction,'G',0);

% objConstant = 0;
[objPoly]=simplifyPolynomial(objPoly1);
if (any(objPoly.supports(1,:),2) == 0 )
    objPoly.coef(1,1) = objConstant + objPoly.coef(1,1);
elseif abs(objConstant) > 1.0e-12
    objPoly.coef = [objConstant; objPoly.coef];
    objPoly.supports = [sparse(1,objPoly.dimVar); objPoly.supports];
    objPoly.noTerms = objPoly.noTerms +1;
end

% objConstant

% ineqPolySys
pointer = 0;
if eqTo2ineqSW == 0
    for i=1:noOfEquations
        pointer = i;
        ineqPolySys{i} = convToPolynomial(noOfVariables,varNames,...
            listOfTerms{i},eqOrIneq{i},rightValue{i});
	%full(ineqPolySys{i}.supports)
	%full(ineqPolySys{i}.coef)
    end
    pointer = noOfEquations;
else % eqTo2ineqSW == 1 
    pointer = 0;
    for i=1:noOfEquations
        if (eqOrIneq{i} == 'G')
            pointer = pointer + 1;
            ineqPolySys{pointer} = convToPolynomial(noOfVariables,varNames,...
                listOfTerms{i},eqOrIneq{i},rightValue{i});
        elseif (eqOrIneq{i} == 'L')
            pointer = pointer + 1;
            ineqPolySys{pointer} = convToPolynomial(noOfVariables,varNames,...
                listOfTerms{i},eqOrIneq{i},rightValue{i});
        else
            pointer = pointer + 1;
            ineqPolySys{pointer} = convToPolynomial(noOfVariables,varNames,...
                listOfTerms{i},'G',rightValue{i});
            pointer = pointer + 1;
            ineqPolySys{pointer} = convToPolynomial(noOfVariables,varNames,...
                listOfTerms{i},'L',rightValue{i}+eqTolerance);
        end
    end
end

% pointer = noOfEquations;
% iEqPolySys --- nonnegativity
if isempty(posVarNames) ~= 1
    ll = size(posVarNames,2);
    for i=1:ll
        lbd = posToInEqPolySys(noOfVariables,varNames,posVarNames{i},lbd);
    end
else
    ll=0;
end

printSW = 0;
if printSW == 1
    nnn = size(lbd,2);
    for i=1:nnn
        fprintf('%+6.2e ',lbd(1,i));
    end
    fprintf('\n');
    for i=1:nnn
        fprintf('%+6.2e ',ubd(1,i));
    end
    fprintf('\n');
end

debug = 0;
if debug
    fprintf('objPoly : ');
    writePolynomials(objPoly);
    writePolynomials(ineqPolySys);
end

fclose(fileIDX);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [statusSW,oneLine] = getOneLine(dataFID)
flowCTRL = 0;
oneLine = '';
while (feof(dataFID)==0) && (flowCTRL== 0)
    inputLine = fgetl(dataFID);
    %	inputLine
    len = length(inputLine);
    if (len > 0) && (inputLine(1)~='*')
        p=1;
        while (p<=len) && (inputLine(p)==' ')
            p = p+1;
        end
        if (p<=len) % & (inputLine(p) ~= '*')
            %			oneLine
            %			inputLine
            % Kojima 11/06/04; to meet MATLAB 5.2
            if isempty(oneLine)
                oneLine = inputLine(p:len);
            else
                oneLine = strcat(oneLine,inputLine(p:len));
            end
            % Kojima 11/06/04; to meet MATLAB 5.2
            %            temp = strfind(inputLine,';');
            temp = findstr(inputLine,';');
            if isempty(temp) == 0
                flowCTRL=1;
            end
        end
    end
end
if flowCTRL==0
    oneLine = '';
    statusSW = -1;
else
    statusSW = 1;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [varNames,p,moreSW] = getListOfNames(oneLine,varNames,p)
while (length(oneLine) > 0)
    [oneName,remLine] = strtok(oneLine,' ,');
    if length(oneName) > 0
        p = p+1;
        varNames{p} = oneName;
    end
    oneLine = remLine;
end
lenLastVar = length(varNames{p});
if varNames{p}(lenLastVar) == ';'
    moreSW = 0;
    if lenLastVar == 1
        p = p-1;
    else
        varNames{p} = varNames{p}(1:lenLastVar-1);
    end
else
    moreSW = 1;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [listOfTerms,eqOrIneq,rightValue] = separate(oneLine)

[formerPart,latterPart] = strtok(oneLine,'=');
eqOrIneq = latterPart(2);
ll = length(latterPart);
rightValue = str2num(latterPart(4:ll-1));
% formerPart
ll = length(formerPart);
k = 0;
while ll > 0
    oneLetter = ' ';
    p = 0;
    while oneLetter == ' '
        p = p+1;
        oneLetter = formerPart(p);
    end
    formerPart = formerPart(p:ll);
    ll = ll - p + 1;
    if ll > 0
        if formerPart(1) == ';'
            ll = 0;
        else
            k = k+1;
            if (formerPart(1) == '-') || (formerPart(1) == '+')
                listOfTerms{k} = formerPart(1);
                formerPart = formerPart(2:ll);
                ll = ll -1;
                p = 0;
                oneLetter = ' ';
                while oneLetter == ' '
                    p = p+1;
                    oneLetter = formerPart(p);
                end
                formerPart = formerPart(p:ll);
                ll = ll - p + 1;
            else
                listOfTerms{k} = '+';
            end
            [oneTerm,formerPart] = strtok(formerPart,'-+;');
            ll = length(formerPart);
            lenOneTerm = length(oneTerm);
            if (lenOneTerm > 0) && ((oneTerm(lenOneTerm) == 'e') || (oneTerm(lenOneTerm) == 'E')) ...
                    && findstr(oneTerm,'.') && (ll > 0) && isempty(findstr(oneTerm,'*'))
                signE = formerPart(1);
                if (signE == '+') || (signE == '-')
                    [oneTerm1,formerPart] = strtok(formerPart,'-+;');
                    oneTerm = strcat(oneTerm,signE);
                    oneTerm = strcat(oneTerm,oneTerm1);
                end
            end
            ll = length(formerPart);
            if isspace(listOfTerms{k})
                listOfTerms{k} = oneTerm;
            else
                listOfTerms{k} = strcat(listOfTerms{k},oneTerm);
            end
        end
    end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [coef,supportVec] = convertOneTerm(noOfVariables,...
    varNames,oneTerm)

supportVec = sparse(1,noOfVariables);
if isempty(oneTerm)
    coef = [];
    supportVec = [];
    return
elseif isletter(oneTerm(2)) 
    if oneTerm(1) == '+'
        coef = 1.0;
    else
        coef = -1.0;
    end
    oneTerm = oneTerm(2:end);
else
    [temp,oneTerm] = strtok(oneTerm,'*');
    coef = str2num(temp);
end
while isempty(oneTerm) ~= 1
    [oneVariable,oneTerm] = strtok(oneTerm,'*');
    kk = length(oneVariable);
    pp = findstr(oneVariable,'^');
    powerPart = 1;
    if isempty(pp) ~= 1
        powerPart = str2num(oneVariable(pp+1:kk));
        oneVariable = oneVariable(1:pp-1);
    end
    i = 1;
    while (i <= noOfVariables)
        if strcmp(oneVariable,varNames{i})
            supportVec(1,i) = supportVec(1,i) + powerPart;
            i = noOfVariables + 1;
        else
            i = i+1;
        end
    end
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function poly = convToPolynomial(noOfVariables,varNames,listOfTerms,...
    eqOrIneq,rightValue)
noOfTerms = size(listOfTerms,2);
if eqOrIneq == 'E'
    poly.typeCone = -1;
else
    poly.typeCone = 1;
end
poly.sizeCone = 1;
poly.degree = 0;
poly.dimVar = noOfVariables;
if abs(rightValue) > 1.0e-10
    poly.noTerms = noOfTerms + 1;
    poly.supports = sparse(1,poly.dimVar);
    poly.coef = -rightValue;
else
    poly.noTerms = noOfTerms;
    poly.supports = [];
    poly.coef = [];
end
for p=1:noOfTerms
    oneTerm = listOfTerms{p};
    [coef,supportVec] = convertOneTerm(noOfVariables,varNames,oneTerm);
    poly.supports = [poly.supports; supportVec];
    poly.coef = [poly.coef;coef];
    %	full(supportVec)
    degree = full(sum(supportVec));
    poly.degree = max(poly.degree,degree);
end
% poly.degree
if eqOrIneq == 'L'
    poly.coef = - poly.coef;
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lbd]  = posToInEqPolySys(noOfVariables,varNames,oneVariable,lbd)

%poly.typeCone = 1;
%poly.sizeCone = 1;
%poly.degree = 1;
%poly.dimVar = noOfVariables;
%poly.noTerms = 1;
%poly.supports = sparse(zeros(1,poly.dimVar));
%poly.coef = [1];

i = 1;
while (i <= noOfVariables)
    if strcmp(oneVariable,varNames{i})
        %		poly.supports(1,i) = 1;
        lbd(i) = max(lbd(i),0);
        i = noOfVariables + 1;
    else
        i = i+1;
    end
end

return

% $Header: /home/waki9/CVS_DB/SparsePOPdev/subPrograms/readGMS.m,v 1.1.1.1 2007/01/11 11:31:50 waki9 Exp $
