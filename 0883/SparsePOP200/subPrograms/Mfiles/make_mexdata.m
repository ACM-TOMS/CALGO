function [ typelist, sizelist, degreelist, dimvarlist, notermslist, supdata, coefdata ] = make_mexdata( objPoly, ineqPolysys )
%[makemexdata]
%This functions convert objPoly and ineqPolysys data
%to special data-type used in mexconv1 function

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
typelist =  objPoly.typeCone ;
sizelist =  objPoly.sizeCone ;
degreelist =  objPoly.degree ;
dimvarlist =  objPoly.dimVar ;
notermslist =  objPoly.noTerms ;

numipoly = size(ineqPolysys,2);

%convert SOCP form to SDP form
for i = 1 : numipoly
    if ineqPolysys{i}.typeCone == 2
        ineqPolysys{i} = SOCPtoSDP(ineqPolysys{i});
    end    
end

%get row , col, val data of  objPoly's coefficients
[sr,sc,sv] = find(objPoly.supports);
if size(objPoly.supports,1) == 1 
    sr = sr';
    sc = sc';
    sv = sv';
end
[cr,cc,cv] = find(objPoly.coef);
if size(objPoly.coef,1 ) == 1
    cr = cr';
    cc = cc';
    cv = cv';
end
if isempty(objPoly.supports)
	totalterms = 0;
	sr = [];
	cr = [];
else	
	totalterms = objPoly.noTerms;
end
%get rows , cols, vals data of  each ineqPolysys's coefficients
for i = 1 : numipoly
    typelist = [ typelist ineqPolysys{i}.typeCone ];
    sizelist = [ sizelist ineqPolysys{i}.sizeCone ];
    degreelist = [ degreelist ineqPolysys{i}.degree ];
    dimvarlist = [ dimvarlist ineqPolysys{i}.dimVar ];
    notermslist = [ notermslist ineqPolysys{i}.noTerms ];
    
    [sr2,sc2,sv2] = find(ineqPolysys{i}.supports);
	if isempty(sr2)
		sr = [];
	end
    [rows,cols] = size(sr2);
    if rows < cols
        sr2 = sr2';
        sc2 = sc2';
        sv2 = sv2';
    end
    sr = [ sr ; sr2 + totalterms];
    sc = [ sc ; sc2 ];
    sv = [ sv ; sv2 ];
    
    [cr2,cc2,cv2] = find(ineqPolysys{i}.coef);
	if isempty(cr2)
		cr = [];
	end
    [rows,cols] = size(cr2);
    if rows < cols
        cr2 = cr2';
        cc2 = cc2';
        cV2 = cV2';
    end
    cr = [ cr ; cr2 + totalterms];
    cc = [ cc ; cc2 ];
    cv = [ cv ; cv2 ];
    
	if ~isempty(ineqPolysys{i}.supports)
    		totalterms = totalterms + ineqPolysys{i}.noTerms;
	end
end

%set coefdata
supdata  = sparse(sr,sc,sv);
coefdata = sparse(cr,cc,cv);

supdata = supdata';
coefdata = coefdata';

% Notice:
% nDim is the number of variables.
% #mono = the number of monomials which appear in the POP
% supdata is nDim times #mono.
%
% If ineqPolySys{1}.typeCone and ineqPolySys{2}.typeCone are 1 or -1
% support of k-th monomial which is in the POP is supdata(:, k),
%
%
% sizeCone = max(objPoly.sizeCone, ineqPolySys{i}.sizeCone)
% #mono = the number of monomials which appear in the POP
% coefdata is #mono \times sizeCone*sizeCone matrix if typeCone > 1 for some ineqPolySys{i}. 
% Otherwise, it is #mono \times sizeCone matrix.
% 
% e.g. ineqPolySys{1}.sizeCone = 2, ineqPolySys{2}.sizeCone = 3 
% If ineqPolySys{1}.typeCone and ineqPolySys{2}.typeCone are 1 or -1
% coef of k-th monomial which is in ineqPolySys{1} is coefdata(1:ineqPolySys{1}.sizeCone, k),
% coef of k-th monomial which is in ineqPolySys{2} is coefdata(:, k),
% where k is a number of monomials in the POP, not in the ineqPolySys{1} or ineqPolySys{2} 
%
% If ineqPolySys{1}.typeCone and ineqPolySys{2}.typeCone are 2 or 3 
% coef of k-th monomial which is in ineqPolySys{1} is coefdata(1:ineqPolySys{1}.sizeCone^2, k),
% coef of k-th monomial which is in ineqPolySys{2} is coefdata(:, k),
% where k is a number of monomials in the POP, not in the ineqPolySys{1} or ineqPolySys{2} 
%
%

return

