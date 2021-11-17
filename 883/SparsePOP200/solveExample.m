function solveExample(probNumbers);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve all POPs in the directories GMSformat and POPformat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% To solve problem 3, 
%   >> solveExample(3);
% To solve problems 2, 8 and 10,  
%   >> solveExample([2, 8, 10]);
% To solve problems 11 through 20,  
%   >> solveExample([11:20]);
% To solve all problems, 
%   >> solveExample;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% List of GAMS scalar format POPs in the directory GMSformat 
% These problems are from GLOBAL Library
%	http://www.gamsworld.org/global/globallib.htm
% Lower and upper bounds of variables are added to some of the problems. 
% Some problems are scaled. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
problemList{1}.name = 'Bex2_1_1.gms';  problemList{1}.relaxOrder = 3; 
problemList{2}.name = 'Bex2_1_2.gms'; 
problemList{3}.name = 'Bex2_1_3.gms'; 
problemList{4}.name = 'Bex2_1_4.gms';  
problemList{5}.name = 'Bex2_1_5.gms';  
problemList{6}.name = 'Bex2_1_8.gms';  
problemList{7}.name = 'Bex3_1_1.gms';  problemList{7}.relaxOrder = 3;
problemList{8}.name = 'Bex3_1_2.gms';  
problemList{9}.name = 'Bex3_1_4.gms'; 
    problemList{9}.relaxOrder = 4; problemList{9}.perturbation = 1.0e-4;
problemList{10}.name = 'Bex5_2_2_case1.gms'; problemList{10}.sparseSW = 0; 
%
problemList{11}.name = 'Bex5_2_2_case2.gms'; problemList{11}.sparseSW = 0; 
problemList{12}.name = 'Bex5_2_2_case3.gms'; problemList{12}.sparseSW = 0; 
problemList{13}.name = 'Bex5_2_5.gms'; problemList{13}.relaxOrder = 1; 
    % nDim = 32, too big to solve.
problemList{14}.name = 'Bex5_3_2.gms'; 
    % nDim = 23, too big to solve.
problemList{15}.name = 'Bex5_4_2.gms'; problemList{15}.relaxOrder = 3; 
problemList{16}.name = 'Bex9_1_1.gms'; problemList{16}.sparseSW = 0; problemList{16}.complementaritySW = 1; 
problemList{17}.name = 'Bex9_1_2.gms'; problemList{17}.sparseSW = 0; problemList{17}.complementaritySW = 1; 
problemList{18}.name = 'Bex9_1_4.gms'; 
    problemList{18}.sparseSW = 0; problemList{18}.perturbation = 1.0e-4; problemList{18}.complementaritySW = 1; 
problemList{19}.name = 'Bex9_1_5.gms'; 
    problemList{19}.sparseSW = 0; problemList{19}.perturbation = 1.0e-4; problemList{19}.complementaritySW = 1; 
problemList{20}.name = 'Bex9_1_8.gms'; problemList{20}.sparseSW = 0; problemList{20}.complementaritySW = 1; 
problemList{21}.name = 'Bex9_2_1.gms'; problemList{21}.sparseSW = 0; problemList{21}.complementaritySW = 1;
%
problemList{22}.name = 'Bex9_2_2.gms'; 
    problemList{22}.sparseSW = 0; problemList{22}.complementaritySW = 1; 
problemList{23}.name = 'Bex9_2_3.gms'; 
    problemList{23}.sparseSW = 0; problemList{23}.complementaritySW = 1; 
problemList{24}.name = 'Bex9_2_4.gms'; 
    problemList{24}.sparseSW = 0; problemList{24}.complementaritySW = 1; 
problemList{25}.name = 'Bex9_2_5.gms'; 
    problemList{25}.sparseSW = 0; problemList{25}.complementaritySW = 1; 
problemList{26}.name = 'Bex9_2_6.gms'; 
    problemList{26}.sparseSW = 0; problemList{26}.complementaritySW = 1; 
problemList{27}.name = 'Bex9_2_7.gms'; 
    problemList{27}.sparseSW = 0; problemList{27}.complementaritySW = 1; 
problemList{28}.name = 'Bex9_2_8.gms'; 
    problemList{28}.sparseSW = 0; problemList{28}.complementaritySW = 1; 
problemList{29}.name = 'Balkyl.gms'; 
    problemList{29}.relaxOrder = 3;
problemList{30}.name = 'Bst_bpaf1a.gms'; 
problemList{31}.name = 'Bst_bpaf1b.gms'; 
%
problemList{32}.name = 'Bst_e05.gms';
problemList{33}.name = 'Bst_e07.gms'; 
problemList{34}.name = 'Bst_jcbpaf2.gms';
problemList{35}.name = 'Bhaverly.gms'; 
problemList{36}.name = 'Babel.gms';
problemList{37}.name = 'alkylation.gms'; 
    problemList{37}.relaxOrder = 3; problemList{37}.perturbation = 1.0e-4;
problemList{38}.name = 'Bst_bpk1.gms';
problemList{39}.name = 'Bst_bpk2.gms'; 
problemList{40}.name = 'Bst_bpv1.gms';
problemList{41}.name = 'Bst_bpv2.gms'; 
%
problemList{42}.name = 'Bst_e33.gms'; 
problemList{43}.name = 'Bst_e42.gms'; 
problemList{44}.name = 'Bst_robot.gms'; 
    problemList{44}.sparseSW = 0; problemList{44}.perturbation = 1.0e-4;
problemList{45}.name = 'meanvar.gms'; problemList{45}.relaxOrder = 1; 
problemList{46}.name = 'mhw4d.gms'; 
problemList{47}.name = 'Bprolog.gms'; 
problemList{48}.name = 'st_cqpjk2.gms';
problemList{49}.name = 'st_e01.gms';
problemList{50}.name = 'st_e09.gms'; problemList{50}.relaxOrder = 3; 
problemList{51}.name = 'st_e10.gms';
%
problemList{52}.name = 'st_e20.gms'; 
problemList{53}.name = 'st_e23.gms'; 
problemList{54}.name = 'st_e24.gms'; 
problemList{55}.name = 'st_e34.gms'; 
problemList{56}.name = 'st_e42.gms'; 
problemList{57}.name = 'st_fp5.gms'; 
problemList{58}.name = 'st_glmp_fp1.gms';
problemList{59}.name = 'st_glmp_fp2.gms'; problemList{59}.relaxOrder = 3; 
problemList{60}.name = 'st_glmp_fp3.gms'; 
problemList{61}.name = 'st_glmp_kk90.gms'; 
%
problemList{62}.name = 'st_glmp_kk92.gms'; 
problemList{63}.name = 'st_glmp_kky.gms'; 
problemList{64}.name = 'st_glmp_ss1.gms';
problemList{65}.name = 'st_glmp_ss2.gms';
problemList{66}.name = 'st_iqpbk1.gms';
problemList{67}.name = 'st_iqpbk2.gms'; 
problemList{68}.name = 'st_jcbpaf2.gms';
problemList{69}.name = 'st_jcbpafex.gms';
problemList{70}.name = 'qp1.gms'; problemList{70}.relaxOrder = 1; 
problemList{71}.name = 'qp2.gms'; problemList{71}.relaxOrder = 1; 
%
problemList{72}.name = 'qp3.gms'; problemList{72}.relaxOrder = 1; 
problemList{73}.name = 'qp5.gms'; problemList{73}.relaxOrder = 1; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% List of the SparsePOP format POPs in the directory POPformat 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
problemList{74}.name = 'Rosenbrock(200,1)'; problemList{74}.printLevel = [1, 1]; 
problemList{75}.name = 'BroydenBand(6)'; problemList{75}.relaxOrder = 3; 
problemList{76}.name = 'BroydenTri(200)'; problemList{76}.printLevel = [1, 1]; 
problemList{77}.name = 'ChainedSingular(200)'; problemList{77}.printLevel = [1, 1]; 
problemList{78}.name = 'ChainedWood(200)'; problemList{78}.printLevel = [1, 1]; 
problemList{79}.name = 'nondquar(200)'; problemList{79}.printLevel = [1, 1]; 
problemList{80}.name = 'nonscomp(8)'; problemList{80}.relaxOrder = 3;
problemList{81}.name = 'optControl(6,2,4,0)'; problemList{81}.printLevel = [1, 1]; 
problemList{82}.name = 'optControl2(200)'; problemList{82}.printLevel = [1, 1]; 
problemList{83}.name = 'randomUnconst(20,2,4,4,3201)'; 
problemList{84}.name = 'randomConst(20,2,4,4,3201)'; problemList{84}.relaxOrder = 3;
problemList{85}.name = 'randomwithEQ(20,2,4,4,3201)'; problemList{85}.relaxOrder = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[param0] = defaultParameter([]);
%param0.mex = 1;
%param0.printFileName = 'result.txt';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 0 
    probNumbers = [1:85]; 
    if param0.symbolicMath == 0
        probNumbers = [1:12, 14:34, 74:85]; 
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


noProblems = length(probNumbers);

for kkk=1:noProblems
    k = probNumbers(kkk); 
	fileName = problemList{k}.name;
    param = param0; 
    if isfield(problemList{k},'sparseSW')
        param.sparseSW = problemList{k}.sparseSW;
    end
    if isfield(problemList{k},'relaxOrder')
        param.relaxOrder = problemList{k}.relaxOrder;
    else
        param.relaxOrder = 2; 
    end
    if isfield(problemList{k},'perturbation')
        param.perturbation = problemList{k}.perturbation;
    end
    if isfield(problemList{k},'complementaritySW')
        param.complementaritySW = problemList{k}.complementaritySW;
    end
    if isfield(problemList{k},'printLevel')
        param.printLevel = problemList{k}.printLevel;
    end
	%{
	s = findstr(fileName,'.gms');
	if isempty(s)
		sdpaName = strcat(fileName,'.dat-s');
	else
		sdpaName = strcat(fileName(1:s),'dat-s');
	end
	param.sdpaDataFile = sdpaName;
	%}
	fprintf('\n\n%d: %s\n', k, fileName);
    [param,SDPobjValue,POP,cpuTime,SeDuMiInfo,SDPinfo] = sparsePOP(fileName,param);
end

return

% $Header: /home/waki9/CVS_DB/SparsePOPdev/example/solveAllGmsProblems.m,v 1.1.1.1 2007/01/11 11:31:49 waki9 Exp $
