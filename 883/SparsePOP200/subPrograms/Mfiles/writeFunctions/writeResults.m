function writeResults(fileId,printLevel,cpuTime,POP,SDPobjValue,SeDuMiInfo,SDPinfo)
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

noOfSDPblocks = length(SDPinfo.SDPblock);
if isfield(SDPinfo,'SOCPblock')
    noOfSOCPblocks = length(SDPinfo.SOCPblock);
else
    noOfSOCPblocks = 0;
end
maxSDPblock = max(SDPinfo.SDPblock);

if noOfSDPblocks > 0
    aveSDPblock = sum(SDPinfo.SDPblock)/noOfSDPblocks;
else
    aveSDPblock = 0;
end

fprintf(fileId,'# Computational results of sparsePOP\n');
fprintf(fileId,'size of A          = [%d, %d]\n',SDPinfo.rowSizeA,SDPinfo.colSizeA);
fprintf(fileId,'no of LP variables = %d\n',SDPinfo.noOfLPvariables);
fprintf(fileId,'no of FR variables = %d\n',SDPinfo.noOfFreevar);
if noOfSOCPblocks > 0
    fprintf(fileId,'no of SOCP blocks  = %d\n',noOfSOCPblocks);
end
fprintf(fileId,'no of SDP blocks   = %d\n',noOfSDPblocks);
fprintf(fileId,'max size SDP block = %d\n',maxSDPblock);
fprintf(fileId,'ave.size SDP block = %6.2e\n',aveSDPblock);

fprintf(fileId,'cpuTime.conversion = %10.2f\n',cpuTime.conversion);
fprintf(fileId,'cpuTime.SeDuMi     = %10.2f\n',cpuTime.SeDuMi);
fprintf(fileId,'cpuTime.total      = %10.2f\n',cpuTime.total);

fprintf(fileId,'SeDuMiInfo.numerr  = %d\n',SeDuMiInfo.numerr);
fprintf(fileId,'SeDuMiInfo.pinf    = %d\n',SeDuMiInfo.pinf);
fprintf(fileId,'SeDuMiInfo.dinf    = %d\n',SeDuMiInfo.dinf);

if SeDuMiInfo.pinf == 0 && SeDuMiInfo.dinf ~= 0
    fprintf(fileId,'Dual SDP problem is infieasible.\n');
elseif SeDuMiInfo.pinf ~= 0 && SeDuMiInfo.dinf == 0
    fprintf(fileId,'Primal SDP problem is infieasible.\n');
    fprintf(fileId,'In this case, the original POP is also infeasible.\n');
elseif SeDuMiInfo.pinf ~= 0 && SeDuMiInfo.dinf ~= 0
    fprintf(fileId,'Primal and Dual SDP problems are infieasible.\n');
    fprintf(fileId,'In this case, the original POP is also infeasible.\n');
else
    fprintf(fileId,'SDPobjValue        = %+17.12e\n',SDPobjValue);
end

if ~isempty(POP.xVect)
    fprintf(fileId,'POP.objValue        = %+17.12e\n',POP.objValue);
    fprintf(fileId,'POP.absError        = %+17.12e\n',POP.absError);
    fprintf(fileId,'POP.scaledError     = %+17.12e\n',POP.scaledError);
end

if printLevel >= 2 && ~isempty(POP.xVect)
    if SeDuMiInfo.pinf == 0 && SeDuMiInfo.dinf == 0
        fprintf(fileId,'POP.xVect = ');
        lenOFx = length(POP.xVect);
        k = 0;
        for j=1:lenOFx
            if mod(k,5) == 0
                fprintf(fileId,'\n');
            end
            k = k+1;
            fprintf(fileId,'%5d:%+17.12e',j,POP.xVect(j));
        end
        fprintf(fileId,'\n');
    end
elseif isempty(POP.xVect)
    if SeDuMiInfo.pinf == 0 && SeDuMiInfo.dinf == 0
        fprintf(fileId,['\n##! SparsePOP can not return approximate\n' ...
            '##! optimal solution information for this problem\n']);
        fprintf(fileId,['##! If you need this solution information,\n' ...
            '##! please solve with param.reduceMomentMatSW = 0\n\n']);
    end
end
return


% $Header: /home/waki9/CVS_DB/SparsePOPdev/subPrograms/writeFunctions/writeResults.m,v 1.1.1.1 2007/01/11 11:31:51 waki9 Exp $
