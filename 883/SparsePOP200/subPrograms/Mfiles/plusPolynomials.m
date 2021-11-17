function polyOut = plusPolynomials(poly1,poly2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% polyOut = poly1 + poly2                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   polynomial.typeCone = -1, 1, 2 0r 3
%   polynomial.sizeCone >= 1
% 	polynomial.dimVar --- a positive integer 
% 	polynomial.degree --- a positive integer 
% 	polynomial.noTerms --- a positive integer 
% 	polynomial.supports --- a noTerms \times nDim matrix
% 	polynomial.coef --- a noTerms \times noTerms vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
if isempty(poly1) && ~isempty(poly2)
    polyOut = poly2;
    return
elseif isempty(poly2) && ~isempty(poly1)    
    polyOut = poly1;
    return    
elseif poly1.dimVar ~= poly2.dimVar
	polyOut = [];
	fprintf('\n!!! poly1.dimVar ~= poly2.dimVar !!!\n\n'); 
	return
elseif poly1.typeCone ~= poly2.typeCone
	polyOut = [];
	fprintf('\n!!! poly1.typeCone ~= poly2.typeCone !!!\n\n'); 
	return
elseif poly1.sizeCone ~= poly2.sizeCone
	polyOut = [];
	fprintf('\n!!! poly1.typeCone ~= poly2.typeCone !!!\n\n'); 
	return
end

polyIn.typeCone = poly1.typeCone;
polyIn.sizeCone = poly1.sizeCone;
polyIn.dimVar = poly1.dimVar;
polyIn.degree = max(poly1.degree,poly2.degree);
polyIn.supports = [poly1.supports; poly2.supports]; 
polyIn.coef = [poly1.coef; poly2.coef]; 
polyIn.noTerms = poly1.noTerms + poly2.noTerms; 
polyOut = simplifyPolynomial(polyIn);
return


% $Header: /home/waki9/CVS_DB/SparsePOPdev/subPrograms/plusPolynomials.m,v 1.1.1.1 2007/01/11 11:31:50 waki9 Exp $
