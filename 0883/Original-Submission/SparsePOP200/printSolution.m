function printSolution(fileId,printLevel,dataFileName,param,...
    SDPobjValue,POP,cpuTime,SeDuMiInfo,SDPinfo)
%
% printSolution
% prints solutions obtained by sparsePOP.
%
% Usage:
%
% printSolution(fileId,printLevel,dataFileName,param,SDPobjValue,...
%	POP,cpuTime,SeDuMiInfo,SDPinfo);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fileId: the file ID where output goes. If this is 1, then the output is
%       the standard output (i.e., the screen). Default is 1.
% printLevel: controls how much information should be printed.
%       A larger value gives more information. Default value is 2.
% dataFileName: the name of the problem.
% The rest of the input arguments must be the outputs of sparsePOP.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% none.

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

if printLevel >= 1
    if (param.SeDuMiSW == 0) | (SDPinfo.infeasibleSW <= -1) | (SDPinfo.infeasibleSW >= 1);
        fprintf(fileId,'\n## Computational Results by sparsePOP.m ##\n');
    else
        fprintf(fileId,'\n## Computational Results by sparsePOP.m');
        fprintf(fileId,' with SeDuMi ##\n');
    end
    fprintf(fileId,'## Printed by printSolution.m ##\n');
end
if ~isempty(dataFileName)
    fprintf(fileId,'# Problem File Name   = %s\n',dataFileName);
end

if printLevel >= 1
    fprintf(fileId,'# parameters:\n');
	fprintf(fileId,'  relaxOrder          = %d\n',param.relaxOrder);
    fprintf(fileId,'  sparseSW            = %d\n',param.sparseSW);
   	if ischar(param.multiCliquesFactor) == 1 
		fprintf(fileId,'  multiCliquesFactror = %s\n',param.multiCliquesFactor);
	elseif isnumeric(param.multiCliquesFactor) == 1
		if param.multiCliquesFactor ~=1 
			fprintf(fileId,'  multiCliquesFactror = %d\n',param.multiCliquesFactor);
		end
	end
	if param.scalingSW == 0
        fprintf(fileId,'  scalingSW           = %d\n',param.scalingSW);
    end
    if param.boundSW == 0 || param.boundSW == 1
        fprintf(fileId,'  boundSW             = %d\n',param.boundSW);
    end
    if param.eqTolerance > 1.0e-10
        fprintf(fileId,'  eqTolerance         = %6.2e\n',param.eqTolerance);
    end
    if param.perturbation > 1.0e-10
        fprintf(fileId,'  perturbation        = %6.2e\n',param.perturbation);
    end
    if param.reduceMomentMatSW == 0
        fprintf(fileId,'  reduceMomentMatSW   = %d\n',param.reduceMomentMatSW);
    end
    if param.complementaritySW == 1 
        fprintf(fileId,'  complementaritySW   = %d\n',param.complementaritySW);
    end
    if param.SeDuMiSW ~= 1
        fprintf(fileId,'  SeDuMiSW            = %d\n',param.SeDuMiSW);
    end
    if ischar(param.SeDuMiOutFile) == 1
        fprintf(fileId,'  SeDuMiOutFile       ');
        fprintf(fileId,'= %s\n',param.SeDuMiOutFile);
    elseif isnumeric(param.SeDuMiOutFile) == 1
		if param.SeDuMiOutFile ~= 0 
			fprintf(fileId,'  SeDuMiOutFile       ');
			fprintf(fileId,'= %d\n',param.SeDuMiOutFile);
		end
    end
    if ~isempty(param.detailedInfFile)
        if ischar(param.detailedInfFile)
            fprintf(fileId,'  detailedInfFile     = %s\n',param.detailedInfFile);
        end
    end
    if ischar(param.sdpaDataFile) && ~isempty(param.sdpaDataFile)
        fprintf(fileId,'  sdpaDataFile        ');
        fprintf(fileId,'= %s\n',param.sdpaDataFile);
    end
    if ischar(param.printFileName) == 1 && ~isempty(param.printFileName)
        fprintf(fileId,'  printFileName       ');
        fprintf(fileId,'= %s\n',param.printFileName);
    elseif isnumeric(param.printFileName) == 1
		if param.printFileName == 0
			fprintf(fileId,'  printFileName       ');
    	    fprintf(fileId,'= %d\n',param.printFileName);
	    end
	end
	%if param.printLevel(1) ~= 2 || param.printLevel(2) ~= 0 
	%		fprintf(fileId,'  printLevel       = [%d, %d]',param.printLevel(1), param.printLevel(2));
	%end
	if param.symbolicMath ~=1 
		fprintf(fileId,'  symbolicMath        = %d\n', param.symbolicMath);
	end
	if param.mex ~= 1 
		fprintf(fileId,'  mex                 = %d\n', param.mex);
	end
end

if (-1 <= SDPinfo.infeasibleSW) && (SDPinfo.infeasibleSW <= 1)
    noOfSDPblocks = length(SDPinfo.SDPblock);
    if noOfSDPblocks > 0
        aveSDPblock = sum(SDPinfo.SDPblock)/noOfSDPblocks;
        maxSDPblock = max(SDPinfo.SDPblock);
    else
        aveSDPblock = 0;
        maxSDPblock = 0;
    end
    if (param.SeDuMiSW == 1) & (SDPinfo.infeasibleSW == 0)
        fprintf(fileId,'# SDP solved by SeDuMi:\n');
    elseif (param.SeDuMiSW == 0) & (SDPinfo.infeasibleSW == 0)
        fprintf(fileId,'# Estimated size of SDP to be solved:\n');
    elseif SDPinfo.reduceAMatSW <= 1 
        fprintf(fileId,'# SDP relaxation problem:\n');
    end
    if SDPinfo.reduceAMatSW <= 1 
        fprintf(fileId,'  size of A           = [%d,%d]\n',SDPinfo.rowSizeA,SDPinfo.colSizeA);
        fprintf(fileId,'  no of nonzeros in A = %d\n',SDPinfo.nonzerosInA);
        fprintf(fileId,'  no of LP variables  = %d\n',SDPinfo.noOfLPvariables);
        fprintf(fileId,'  no of FR variables  = %d\n',SDPinfo.noOfFreevar);
        if SDPinfo.SOCPblock > 0
            fprintf(fileId,'  no of SOCP blocks    = %d\n',SDPinfo.SOCPblock);
        end
        fprintf(fileId,'  no of SDP blocks    = %d\n',noOfSDPblocks);
        fprintf(fileId,'  max size SDP block  = %d\n',maxSDPblock);
        fprintf(fileId,'  ave.size SDP block  = %6.2e\n',aveSDPblock);
    end
end
if ((param.SeDuMiSW == 1) && (SDPinfo.infeasibleSW == 0)) | (SDPinfo.infeasibleSW == -1)
    if ~isempty(SeDuMiInfo)
        fprintf(fileId,'# SeDuMi information:\n');
        fprintf(fileId,'  SeDuMi.pars.eps     = %6.2e\n',param.SeDuMiEpsilon);
        fprintf(fileId,'  SeDuMiInfo.numerr   = %d\n',SeDuMiInfo.numerr);
        if SeDuMiInfo.numerr == 1
            fprintf(fileId,'      SeDuMi stopped before SeDuMi.pars.eps = %6.2e is attained.\n',param.SeDuMiEpsilon); 
        elseif SeDuMiInfo.numerr >= 2
            fprintf(fileId,'      SeDuMi stopped because of serious numerical difficulties.\n'); 
        end
        fprintf(fileId,'  SeDuMiInfo.pinf     = %d\n',SeDuMiInfo.pinf);
        if SeDuMiInfo.pinf ~= 0 
            fprintf(fileId,'      Primal SDP is infeasible.\n');
            fprintf(fileId,'      No finite lower bound has been computed for POP objective function.\n');
        end 
        fprintf(fileId,'  SeDuMiInfo.dinf     = %d\n',SeDuMiInfo.dinf);
        if SeDuMiInfo.dinf ~= 0
            fprintf(fileId,'      Dual SDP problem is infeasible.\n');
            fprintf(fileId,'      POP is probably infeasible.\n');
        end
    end
    if (SDPinfo.infeasibleSW == -1) | (SeDuMiInfo.pinf == 0 && SeDuMiInfo.dinf == 0)
        fprintf(fileId,'# Approximate optimal value information:\n');
        if SDPinfo.infeasibleSW == 0 
            fprintf(fileId,'  SDPobjValue         = %+13.7e\n',SDPobjValue);
        end
        if ~isempty(POP.xVect)
            if issparse(POP.objValue)
                POP.objValue = full(POP.objValue);
            end
			fprintf(fileId,'  POP.objValue        = %+13.7e\n',POP.objValue);
            relobj = abs(POP.objValue-SDPobjValue)/max(1,abs(POP.objValue));
			if SDPinfo.infeasibleSW == 0 
                fprintf(fileId,'  relative obj error  = %+8.3e\n',relobj);
            end
            if issparse(POP.absError)
                POP.absError = full(POP.absError);
            end
            fprintf(fileId,'  POP.absError        = %+8.3e\n',POP.absError);
            if issparse(POP.scaledError)
                POP.scaledError = full(POP.scaledError);
            end
            fprintf(fileId,'  POP.scaledError     ');
            fprintf(fileId,'= %+8.3e\n',POP.scaledError);
        end
%     else
%         if SeDuMiInfo.pinf == 0 && SeDuMiInfo.dinf ~= 0
%            fprintf(fileId,'  Dual SDP problem is infieasible.\n');
%        elseif SeDuMiInfo.pinf ~= 0 && SeDuMiInfo.dinf == 0
%            fprintf(fileId,'  Primal SDP problem is infieasible.\n');
%			fprintf(fileId,'  Increase param.relaxOrder and solve this POP.\n');
%        elseif SeDuMiInfo.pinf ~= 0 && SeDuMiInfo.dinf ~= 0
%            fprintf(fileId,'  Primal and Dual SDP problems are infieasible.\n');
%            fprintf(fileId,'  The original POP may be also infeasible.\n');
%        end
    end
elseif (SDPinfo.infeasibleSW == 1)
    fprintf(fileId,'# Primal SDP is infeasible!\n');
    fprintf(fileId,'  No finite lower bound has been computed for the POP objective function!\n');
elseif (SDPinfo.infeasibleSW == 2)
    fprintf(fileId,'# POP is infeasible!\n');
elseif (SDPinfo.infeasibleSW == -2)
    fprintf(fileId,'# Approximate optimal value information:\n');
    fprintf(fileId,'  POP has a unique feasible solution\n');
	if issparse(POP.objValue)
		POP.objValue = full(POP.objValue);
	end
    fprintf(fileId,'  POP.objValue        = %+13.7e\n',POP.objValue);
	if issparse(POP.absError)
		POP.absError = full(POP.absError);
	end
    fprintf(fileId,'  POP.absError        = %+8.3e\n',POP.absError);
    if issparse(POP.scaledError)
        POP.scaledError = full(POP.scaledError);
    end
    fprintf(fileId,'  POP.scaledError     ');
    fprintf(fileId,'= %+8.3e\n',POP.scaledError);
end

if (param.SeDuMiSW == 1) & (SDPinfo.infeasibleSW == 0)
    fprintf(fileId,'# cpu time:\n');
    fprintf(fileId,'  cpuTime.conversion  = %8.2f\n',cpuTime.conversion);
    fprintf(fileId,'  cpuTime.SeDuMi      = %8.2f\n',cpuTime.SeDuMi);
    fprintf(fileId,'  cpuTime.total       = %8.2f\n',cpuTime.total);
else % param.SeDuMiSW == 0
    fprintf(fileId,'# cpu time:\n');
    fprintf(fileId,'  cpuTime.conversion  = %8.2f\n',cpuTime.conversion);
    fprintf(fileId,'  cpuTime.total       = %8.2f\n',cpuTime.conversion);    
end

if (SDPinfo.infeasibleSW <= 0) && (param.SeDuMiSW==1) ...
   && (printLevel >= 2) && ~isempty(POP.xVect) % && ~isempty(SeDuMiInfo) 
    fprintf(fileId,'# Approximate optimal solution information:\n');
    fprintf(fileId,'  POP.xVect = ');
    lenOFx = length(POP.xVect);
    k = 0;
    for j=1:lenOFx
        if mod(k,5) == 0
            fprintf(fileId,'\n  ');
        end
        k = k+1;
        fprintf(fileId,'%4d:%+13.7e ',j,POP.xVect(j));
    end
    fprintf(fileId,'\n');
elseif (param.SeDuMiSW==1) && (printLevel >= 2) && isempty(POP.xVect) && ~isempty(SeDuMiInfo)
    if SeDuMiInfo.pinf == 0 && SeDuMiInfo.dinf == 0
        fprintf(fileId,['\n##! SparsePOP can not return approximate\n' ...
            '##! optimal solution information for this problem\n']);
        fprintf(fileId,['##! If you need this solution information,\n' ...
            '##! please solve with param.reduceMomentMatSW = 0\n\n']);
    end
end
return

% $Header: /home/waki9/CVS_DB/SparsePOPdev/subPrograms/writeFunctions/printSolution.m,v 1.4 2007/01/16 09:11:40 waki9 Exp $
