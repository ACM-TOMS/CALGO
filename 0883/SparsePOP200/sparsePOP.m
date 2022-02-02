function [param,SDPobjValue,POP,cpuTime,SeDuMiInfo,SDPinfo] = ... 
    sparsePOP(objPoly,ineqPolySys,lbd,ubd,param)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SPARSEPOP  a SPARSE SOS and SDP relaxations to a POP.
% 
%  GENERAL DESCRIPTION
% 
% 	POP is an abbreviation of Polynomial Optimization Problem.
% 	Calling this function, one obtains the optimal value of an SDP relaxation
% 	problem for a POP. A typical invoking line may be:
% 
% >>  [param,SDPobjValue,POP,cpuTime,SeDuMiInfo,SDPinfo] = sparsePOP(DataFile);
% 
% 	The meanings of argument and return values are described below.
% 
%  FILE ARGUMENT
% 
% 	DataFile must be a string containing the file name. The DataFile must
% 	be written in the GMS format or POP format. In the case of GMS format,
% 	for example, call with
% >> sparsePOP('Bex3_1_1.gms');
% 	Don't forget the extension .gms; the file name must be exact.
% 	In the case of POP format, you should write it as if you are calling matlab
% 	function, e.g.,
% >> sparsePOP('BroydenTri(10)');
% 	See userGuide.pdf for more details.
% 
%  RETURN VALUES
% 
% 	param is a structure of parameters used in the execution.
% 	For more details on param structure, see below.	
% 
% 	SDPobjValue is the optimal value of the SDP relaxation problem.
% 
% 	POP is a structure containing information on the POP. Specifically,
% 	POP.xVect is a tentative solution for POP  calculated by the SDP relaxation.
% 	POP.objvalue is the objective value of POP.xVect.
% 	POP.absError is the maximum feasibility violation of POP.xVect.
% 	Their scaled values are also stored in POP.objValScaled and POP.scaledError.
% 
% 	cpuTime is the time consumed by the program execution.
% 	cpuTime.SeDuMi is the time consumed by SeDuMi (SDP Solver),
% 	cpuTime.conversion is the time needed for generating SDP relaxation from POP,
% 	and cpuTime.total is the total.
% 	
% 	SeDuMiInfo is the information passed by SeDuMi. See the manual of SeDuMi
% 	for its details.
% 
% 	SDPinfo contains some statistics of the SDP relaxation problem.
% 
%  OPTIONAL PARAMETERS
% 
% 	You can pass additional parameters in the second argument.
% 
% >>  [param,SDPobjValue,POP,cpuTime,SeDuMiInfo,SDPinfo] = ...
%   sparsePOP(DataFile);
% 
% 	Below are some of the entries of param frequently used. See userGuide.pdf
% 	for the complete information.
% 
% 	param.relaxOrder	Relaxation order for the SDP relaxation.
% 
% 	param.sparseSW		if 1, sparse SDP relaxation is used.
% 		if 0, then the dense (Lasserre's original) SDP relaxation is used.
% 
% 	param.perturbation 	If 1, then sparsePOP perturbs the objective function
% 		so that the resulting problem has a unique optimal solution.
% 		If 0, no perturbation is performed.
% 
% 	param.symbolicMath	if 1,   symbolic math toolbox is used.
% 		Be sure that you have	purchased the symbolic math toolbox
% 		from the MathWorks. if 0, it is not used.
% 
%  ANOTHER ARGUMENT STYLE
% 
% 	Invoking sparsePOP by the following line:
% 
% >> [param,SDPobjValue,POP,cpuTime,SeDuMiInfo,SDPinfo] = ...
%     sparsePOP(objPoly,ineqPolySys,lbd,ubd,param)
% 
% 	one can directly pass all the information of POP through MATLAB structures.
% 	See userGuide.pdf for the description of each component of the arguments.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

fprintf('\nSparsePOP 2.00 by H.Waki, S.Kim, M.Kojima,');
fprintf(' M.Muramatsu and H.Sugimoto\n');
fprintf('                                            ');
fprintf('                June 2007\n\n');
% Check whether SeDuMi has been already installed or not.
if exist('sedumi') ~= 2
    error('User must install SeDuMi and add it to your MATLAB search path before you use SparsePOP.');
end

% Check whether the input problemData is a gms file or
% polynomial format file and set problemName by eleminating '.gms' from
% the gms file or '(...)' from the polynomial format file

gmsSW = 0;
polySW = 0;
mFileSW = 0;
if nargin == 1
    problemData = objPoly;
    param = [];
elseif nargin == 2
    problemData = objPoly;
    param = ineqPolySys;
elseif nargin == 3
    error('Input has something wrong.');
else
    polySW = 1;
    if nargin == 4
    	param = [];
    end
	if isfield(objPoly,'dimVar')
		nDim = objPoly.dimVar;
	else
		error('Set a value in the field dimVar of objPoly.');
	end
	mDim = size(ineqPolySys,2);
	nDim = num2str(nDim);
	mDim = num2str(mDim);
	problemData = strcat('nDim: ', nDim, ', mDim: ', mDim, '.'); 
end

if polySW == 0
    % Input is described in either the GAMS format or the SparsePOP format.
    gmsForm = strfind(problemData,'.gms');
    if length(gmsForm) == 1
        %%
        %% Input is a gms file.
        %%
        gmsSW = 1;
    elseif length(gmsForm) >1 || isempty(problemData)
        %%
        %% if input has more than one string 'gms', we regard it as error.
        %%
        error('Input problem must be a gms file or polynomial format file.');
    elseif isempty(gmsForm) && ~isempty(problemData)
        %%
        %% Input is an m-file which returns POP in the SparsePOP format.
        %%
        mFileSW = 1;
    end
end


% param
param = defaultParameter(param);

% read Data
if gmsSW == 1 % the input file is a gms file
    [objPoly,ineqPolySys,lbd,ubd] = readGMS(problemData,param.symbolicMath);
elseif mFileSW == 1
    [objPoly,ineqPolySys,lbd,ubd] = eval(problemData);
end


% Add objPoly.sizeCone = 1 if it is not specified by the user
if ~isfield(objPoly,'sizeCone')
    objPoly.sizeCone = 1;
end
% Add ineqPolySys{i}.sizeCone = 1 if it is not specified by the user
for i=1:size(ineqPolySys,2);
    if ~isfield(ineqPolySys{i},'sizeCone')
         ineqPolySys{i}.sizeCone = 1;
    end
end

if isfield(param,'multiCliquesFactor') && ischar(param.multiCliquesFactor) 
    param.multiCliquesFactor = objPoly.dimVar;
end

% check inputs
if issparse(lbd)
	lbd = full(lbd);
end
if issparse(ubd)
	ubd = full(ubd);
end

continueSW = checkPOP(objPoly,ineqPolySys,lbd,ubd,param);
if continueSW == 0
    error('## Some inconsistensy in input data, objPoly, ineqPolySys, lbd, ubd and param ##\n'); 
end

% Compute relaxOrder
rOtmp = ceil(objPoly.degree/2);
for i=1:size(ineqPolySys,2)
	tmpdeg = ceil(ineqPolySys{i}.degree/2);
	rOtmp = max(rOtmp,tmpdeg);
end
if ~isfield(param,'relaxOrder')
    param.relaxOrder = rOtemp;
elseif isfield(param,'relaxOrder') && param.relaxOrder < rOtmp 
	param.relaxOrder = rOtmp;
end

if param.mex == 1
    %%
    %% If mex files exist in your MATLAB search path, you uese mex version  
    %%
    [param,SDPobjValue,POP,cpuTime,SeDuMiInfo,SDPinfo] = ...
        SDPrelaxationMex(param,objPoly,ineqPolySys,lbd,ubd);
elseif param.mex == 0
    %%
    %% If mex files do not exist in your MATLAB search path, you uese no
    %% mex version
    %%
    [param,SDPobjValue,POP,cpuTime,SeDuMiInfo,SDPinfo] = ...
        SDPrelaxation(param,objPoly,ineqPolySys,lbd,ubd);
else
	error('set param.mex = 1 or 0');
end
if param.printLevel(1) >= 1
    printSolution(1,param.printLevel(1),problemData,param,SDPobjValue,...
        POP,cpuTime,SeDuMiInfo,SDPinfo);
end

if ischar(param.printFileName) && ~isempty(param.printFileName)
	if param.printLevel(2) == 0
		param.printLevel(2) = 2;
	end 
    fileId = fopen(param.printFileName,'a+');  
    printSolution(fileId,param.printLevel(2),problemData,param,...
        SDPobjValue,POP,cpuTime,SeDuMiInfo,SDPinfo); 
    fclose(fileId); 
end


return
