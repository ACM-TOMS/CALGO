function param = defaultParameter(param)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting default parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Details of Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%
% 1. Parameters to control the basic relaxation scheme.
%
% param.relaxOrder = ¥omega = the relaxation order;  
%   Default value = the minimum relaxation order ¥omega_{¥max}.
%
% param.sparseSW = 1 if you use sparse relaxation;  
%                  0 for dense relaxation; 
%	Default value = 1.
%
% param.multiCliqueFactor 
%   = 0 for no expansion of cliques; 
%   = 1 for combining cliques as long as their sizes do not exceed 
%       the maximum size of all maximal cliquees; 
%   = 'objPoly.dimVar' for combining cliques as long as possible;  
%	Default value = 1.
%
%%%%%%%%%%
% 2. Switch to handle numerical difficulties. 
%
% param.scalingSW
%	= 0 for no scaling.  
%	= 1 for scaling; 
% 	Default value = 2; 
%
% param.boundSW
%   = 0 for no bounds for any y_{\alpha}; 
%   = 1 for bounds for all y_{\alpha};
%   = 2 for bounds for all y_{\alpha} and eliminating redundant bounds; 
% 	Default value = 2; 
%
% param.eqTolerance
%       Convert one equality into two inequality; 
%       if 1.0e-10 \leq param.eqTolerance then f(x) = 0 ===>  f(x) \geq 0  and
%        -f(x) \geq -param.eqTolerance. 
%       Else keep one equality as it is. 
%   Default value = 0.0; 
%
% param.perturbation 
%   = 0 for no pertubation to the objective polynomial; 
%   = 1.0e-5 for a perturbation to the objective polynomial with p, 
%     |p_i| <= 1.0e-5.
% 	Default value = 0; 
%
% param.reduceMomentMatSW 
% 	= 0 for no reduction of moment matrices; 
%   = 1 for reduction of moment matrices by eliminating 
%       redundant elements in moment matrices; 
%   Default value= 1.
%
% param.complementaritySW 
%       If x_ixj = 0 is invloved in equality constraits, 
%       any variable y_{\alpha} correspoinding to a monomial x^{\alpha} 
%       such that \alpha_i \geq 1 and \alpha_j \geq 1 is set to be zero 
%       and eliminaed from the relaxed problem. 
%       Set param.complementaritySW = 1 only when the complementarity condition 
%       is involved in constraints. 
%   = 0 for no reduction in moment matrices using complementarity;
%   = 1 for reduction in moment matrices using complementarity;
%   Default value = 0; 
% 
% param.reduceAMatSW
%       If param.reduceAMatSW = 1, then
%       (a) Eliminate some fixed variables from a POP before applying the 
%           sparse SDP relaxation,  
%       (b) When the equality constraints of the SeDuMi format primal SDP
%           are linearly dependent, eliminated some equalities to restore 
%           the linear independence.
%   = 0 for no (a) and (b);
%   = 1 for (a) and (b);
% Defaulat value = 1;
%
%%%%%%%%%%
% 3. Parameters for SDP solvers.
%
% param.SeDuMiSW 
% 	= 1 for solving SDP by SeDuMi;  
% 	= 0 only for displaying information on the SDP to be solved; 
%   Default value = 1.
%
% param.SeDuMiEpsilon
%       A stopping criteria for the duality gap in SeDuMi;
%   = any nonnegative real number;
%   Default value = 1.0e-9. 
%
% param.SeDuMiOutFile 
%       Specifies where SeDuMi output goes.
% 	= 1 for the standard output (screen)
% 	= 0 for no output
%   = 'filename' for output file name
% Default value = [].
%
% param.sdpaDataFile
%       Specify SDPA sparse format data such that param.sdpaDataFile =
%       'fileName.dat-s', for example, param.sdpaDataFile = 'test.dat-s'; 
%   = [] for no SDPA sparse format data output; 
%   = 'fileName.dat-s';
% 	Default value = []; 
% 
%%%%%%%%%%%
% 4. Parameters for printing numerical results. 
%
% param.detailedInfFile 
%   = 0 for no output of detailed information;
%   = 1 for the screen output of detailed information;
%   = 'filename' for output file name of detailed information;
%   Default value = [];
%
% param.printFileName 
%   = 0 for no solution information; 
%   = 1 for solution informatin;
%   = 'filename' for output file name;
%   Default = 1.
%
% param.printLevel = [a,b], 
%   where a is for display out put, and b is for the file output. 
%   a = 0 for no information on the computational result. 
%       1 for some informtion without an optimal solution. 
%       2 for detailed solution information. 
%   b = 0 for no information on the computational result. 
%       1 for some informtion without an optimal solution. 
%       2 for detailed solution information. 
%   Default = [2, 2].
%
%%%%%%%%%%%
% 5. Parameters to use Symbolic Math Toolbox and C++ subrouties
%
% param.symbolicMath
%   = 1 to use Symbolic MathToolbox;
%   = 0 otherwise;
%   Default value 
%       = 1 if Symbolic Math Toolbox is available;
%       = 0 otherwise; 
% 
% param.mex
%	= 1 to use C++ subrouties;
%	= 0 otherwise;
%	Defaule value
%	if exist('mexconv1') == 3 && exist('mexconv2') == 3 
%    	param.mex = 1;
%	else
%		param.mex = 0;
%	end
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

%%%%%%%%%%%
% 1. Parameters that control the basic relaxation scheme.
if ~isfield(param,'relaxOrder')
    param.relaxOrder = 1; 
%   param.relaxOrder will be updated to
%       max{the minimum relaxation order ¥omega_{¥max}, param.relaxOrder}.
end

if ~isfield(param,'sparseSW')
    param.sparseSW = 1; 
%   Default:
%   param.sparseSW = 1; 
end

if ~isfield(param,'multiCliquesFactor')
    param.multiCliquesFactor = 1; 
%   Default:
%   param.multiCliquesFactor = 1; 
end

%%%%%%%%%%
% 2. Switch to handle numerical difficulties. 
%
if ~isfield(param,'scalingSW') 
	param.scalingSW = 1;
%   Default:
%	param.scalingSW = 1; 
end

if ~isfield(param,'boundSW') 
	param.boundSW = 2; 
%   Default:
%	param.boundSW = 2; 
end

if ~isfield(param,'eqTolerance') 
	param.eqTolerance = 0.0; 
end

if ~isfield(param,'perturbation')
    param.perturbation = 0.0; 
%   Default
%   param.perturbation = 0.0; 
end

if ~isfield(param,'reduceMomentMatSW') 
	param.reduceMomentMatSW = 1;
%   Default:
%	param.reduceMomentMatSW = 1; 
end

if ~isfield(param,'complementaritySW') 
	param.complementaritySW = 0;
%   Default:
%	param.complementaritySW = 0;
end

if ~isfield(param,'reduceAMatSW') 
	param.reduceAMatSW= 1;
%   Default:
%	param.reduceAMatSW = 1;
end

%%%%%%%%%%
% 3. Parameters for SDP solvers
%

if ~isfield(param,'SeDuMiSW')
    param.SeDuMiSW = 1; 
%   Default:
%   param.SeDuMiSW = 1; 
end

if ~isfield(param,'SeDuMiEpsilon')
    param.SeDuMiEpsilon = 1.0e-9;
%   Default:
%   param.SeDuMiEpsilon = 1.0e-9;
end

if ~isfield(param,'SeDuMiOutFile')
    param.SeDuMiOutFile = 0;
%   Default:
%   param.SeDuMiOutFile = 0;
end

if ~isfield(param,'sdpaDataFile') 
	param.sdpaDataFile = '';
%   Default:
%	param.sdpaDataFile = '';
end

%%%%%%%%%%%
% 4. Parameters for printing numerical results. 

if ~isfield(param,'detailedInfFile')
	param.detailedInfFile = '';
%   Default:
%	param.detailedInfFile = '';
end

if ~isfield(param,'printFileName')
    param.printFileName = 1; 
%   Default:
%   param.printFileName = 1; 
end

if ~isfield(param,'printLevel')
    param.printLevel = [2, 2]; 
    if param.printFileName == 0
       param.printLevel(2) = 0; 
    elseif param.printFileName == 1
       param.printLevel(2) = 0; 
    end
%   Default:
%   param.printLevel = [2, 2]; 
end
    
%%%%%%%%%%%
% 5. Parameters to use Symbolic Math Toolbox and C++ subrouties

if ~isfield(param,'symbolicMath')
%   Default:
	A = ver('Symbolic');
	if ~isempty(A)
		x = strfind(A.Name, 'Symbolic Math Toolbox');
		if ~isempty(x)
			param.symbolicMath = 1;
		else
			param.symbolicMath = 0;
		end
	else
			param.symbolicMath = 0;
	end
end

if ~isfield(param,'mex')
	if exist('mexconv1') == 3 && exist('mexconv2') == 3 
    	param.mex = 1;
	else
		param.mex = 0;
	end
elseif param.mex == 1
	if exist('mexconv1') ~= 3 || exist('mexconv2') ~= 3
    	%fprintf('*** C++ files mexconv1.cpp or mexconv2.cpp has not been compiled by compileSparsePOP.m. ***\n');
    	fprintf('*** mexconv1 or mexconv2 compiled by comileSparsePOP.m can not be found. ***\n');
    	fprintf('*** SparsePOP sets param.mex = 0. *** \n');
		param.mex = 0;
	end
end

return
