function writeParameters(fileId,param)
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

fprintf(fileId,'# parameters:\n');
fprintf(fileId,'  relaxOrder         = %d\n',param.relaxOrder);
fprintf(fileId,'  sparseSW           = %d\n',param.sparseSW);
fprintf(fileId,'  multiCliquesFactor = %6.2e\n',param.multiCliquesFactor);
fprintf(fileId,'  scalingSW          = %d\n',param.scalingSW);
fprintf(fileId,'  boundSW            = %d\n',param.boundSW);
fprintf(fileId,'  eqTolerance        = %6.2e\n',param.eqTolerance);
fprintf(fileId,'  perturbation       = %6.2e\n',param.perturbation);
fprintf(fileId,'  reduceMomentMatSW  = %d\n',param.reduceMomentMatSW);
fprintf(fileId,'  complementaritySW  = %d\n',param.complementaritySW);
fprintf(fileId,'  SeDuMiSW           = %d\n',param.SeDuMiSW);
if ischar(param.SeDuMiOutFile)
        fprintf(fileId,'  SeDuMiOutFile      = %s\n',param.SeDuMiOutFile);
elseif isnumeric(param.SeDuMiOutFile)
        fprintf(fileId,'  SeDuMiOutFile      = %d\n',param.SeDuMiOutFile);
end
if ischar(param.detailedInfFile)
        fprintf(fileId,'  detailedInfFile    = %s\n',param.detailedInfFile);
elseif isnumeric(param.detailedInfFile)
        fprintf(fileId,'  detailedInfFile    = %d\n',param.detailedInfFile);
end
if ischar(param.sdpaDataFile)
        fprintf(fileId,'  sdpaDataFile       = %s\n',param.sdpaDataFile);
else
        fprintf(fileId,'  sdpaDataFile       = ''\n');
end
if ischar(param.printFileName)
        fprintf(fileId,'  printFileName      = %s\n',param.printFileName);
elseif isnumeric(param.printFileName)
        fprintf(fileId,'  printFileName      = %d\n',param.printFileName);
end
fprintf(fileId,'  printLevel         = [%d, %d]\n',param.printLevel(1),param.printLevel(2));
fprintf(fileId,'  SeDuMiEpsilon      = %6.2e\n',param.SeDuMiEpsilon);
fprintf(fileId,'  symbolicMath       = %d\n',param.symbolicMath);
fprintf(fileId,'  mex                = %d\n\n',param.mex);
return


% $Header: /home/waki9/CVS_DB/SparsePOPdev/subPrograms/writeFunctions/writeParameters.m,v 1.1.1.1 2007/01/11 11:31:51 waki9 Exp $
