function [objPoly,ineqPolySys,lbd,ubd] = example1;
%%%%%%%%%%%%%%
% example1.m  
%%%%%%%%%%%%%%
%
% The SparsePOP format data for the example1: 
%
% minimize -2*x1 +3*x2 -2*x3
% subject to 
%       x1^2 + 3*x2^2 -2*x2*x3 +3*x3^2 -17*x1 +8*x2 -14*x3 >= -19, 
%       x1 + 2*x2 + x3 <= 5, 
%       5*x2 + 2*x3 = 7, 
%       0 <= x1 <= 2, 0 <= x2 <= 1. 
% 
% To solve the problem by sparsePOP.m: 
% >> param.relaxOrder = 3;
% >> sparsePOP('example1',param);
% 
% This problem is also described in terms of the GAMS scalar format in the 
% file example1.gms.  See Section 3 of the manual.
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

%'example1'
% objPoly 
% -2*x1 +3*x2 -2*x3
        objPoly.typeCone = 1;
        objPoly.dimVar   = 3;    
        objPoly.degree   = 1;    
        objPoly.noTerms  = 3; 
        objPoly.supports = [1,0,0; 0,1,0; 0,0,1]; 
        objPoly.coef     = [-2; 3; -2]; 
% ineqPolySys
% 19 -17*x1 +8*x2 -14*x3 +6*x1^2 +3*x2^2 -2*x2*x3 +3*x3^2 >= 0,
        ineqPolySys{1}.typeCone = 1;
        ineqPolySys{1}.dimVar   = 3;    
        ineqPolySys{1}.degree   = 2;    
        ineqPolySys{1}.noTerms  = 8; 
        ineqPolySys{1}.supports = [0,0,0; 1,0,0; 0,1,0; 0,0,1; ...
                                   2,0,0; 0,2,0; 0,1,1; 0,0,2]; 
        ineqPolySys{1}.coef     = [19; -17; 8; -14; 6; 3; -2; 3]; 
%    
% 5 -x1 -2*x2 -x3  >= 0.   
        ineqPolySys{2}.typeCone = 1;
        ineqPolySys{2}.dimVar   = 3;    
        ineqPolySys{2}.degree   = 1;    
        ineqPolySys{2}.noTerms  = 4; 
        ineqPolySys{2}.supports = [0,0,0; 1,0,0; 0,1,0; 0,0,1]; 
        ineqPolySys{2}.coef     = [5; -1; -2; -1];
%    
% 7 -5*x2 -2*x3 = 0. 
        ineqPolySys{3}.typeCone = -1;
        ineqPolySys{3}.dimVar   = 3;
        ineqPolySys{3}.degree   = 1;
        ineqPolySys{3}.noTerms  = 3;
        ineqPolySys{3}.supports = [0,0,0; 0,1,0; 0,0,1];
        ineqPolySys{3}.coef     = [7; -5; -2];
% lower bounds for variables x1, x2 and x3. 
% 0 <= x1, 0 <= x2, -infinity < x3: 
        lbd = [0,0,-1.0e10];
% upper bounds for variables x1, x2 and x3
% x1 <= 2, x2 <= 1, x3 < infinity: 
        ubd = [2,1,1.0e10];
return
% end of example1.m
