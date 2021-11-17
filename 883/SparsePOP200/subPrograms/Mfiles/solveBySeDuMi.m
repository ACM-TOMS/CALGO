function [x,y,SeDuMiInfo] = solveBySeDuMi(fileId,A,b,c,K,param)

% conversion from the SDPA sparse format data into the SeDuMi
% format data
% Primal:
%	minimize c^T x
%	sub. to  A x = b, x \succeq 0.
% Dual:
%	maximize b^T y
%	sub. to  c - A^T y \succeq 0
% PSDP is converted into the dual problem.


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
if fileId > 0
	printLevel = 0;
    writeSeDuMiInputData(fileId,printLevel,A,b,c,K);
end

if isfield(param,'SeDuMiOutFile')
    if ischar(param.SeDuMiOutFile) && ~isempty(param.SeDuMiOutFile)
        %% Output to file
        pars.fid = fopen(param.SeDuMiOutFile,'a+');
        fprintf(pars.fid,'\n');
    elseif isnumeric(param.SeDuMiOutFile) && param.SeDuMiOutFile > 0
        %% Output to screen
        pars.fid = 1;
    else
        %pars.fid = 0;
        pars.fid = '';
    end
else
    %% No output
    %pars.fid = 0;
    pars.fid = '';
end
% applying the SeDuMi to the SDP
pars.eps = param.SeDuMiEpsilon;
%save mat
%pars.errors = 1;
[x,y,SeDuMiInfo] = sedumi(A,b,c,K,pars);
if isfield(param,'SeDuMiOutFile') && ischar(param.SeDuMiOutFile) && ~isempty(param.SeDuMiOutFile)
    fclose(pars.fid);
end

return
