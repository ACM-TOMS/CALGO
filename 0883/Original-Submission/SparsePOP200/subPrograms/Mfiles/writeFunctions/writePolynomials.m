function writePolynomials(fileId,polynomials)
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
noOfPolynomials = size(polynomials,2); 
% if 1 == noOfPolynomials 
% Kojima 2006/04/06
% To handle an cell array of a single polynomial
if ~iscell(polynomials)
	%	fprintf(fileId,'\n***** Polynomial *****\n\n'); 
	fprintf(fileId,'typeCone = %+d ',polynomials.typeCone); 
	fprintf(fileId,'sizeCone = %d ',polynomials.sizeCone); 
	fprintf(fileId,'dimVar = %d  ',polynomials.dimVar); 
	fprintf(fileId,'degree = %d  ',polynomials.degree); 
	fprintf(fileId,'noTerms = %d\n',polynomials.noTerms); 
	fprintf(fileId,'supportSet = \n'); 
	for i=1:polynomials.noTerms
		fprintf(fileId,'%3d:',i); 
		for j=1:polynomials.dimVar
			fprintf(fileId,' %2d',polynomials.supports(i,j));
		end
		fprintf(fileId,'\n');
	end
	fprintf(fileId,'coefficient = \n'); 
	if polynomials.sizeCone == 1
		k = 0;
		for i=1:polynomials.noTerms
			k = k+1; 
			if mod(k,10) == 0
				fprintf(fileId,'\n');
				k = 0;
			end
			fprintf(fileId,'%3d:%+6.2e ',i,polynomials.coef(i,1));
		end
    	fprintf(fileId,'\n');
	else
		colSize = size(polynomials.coef,2); 
		for i=1:polynomials.noTerms
			fprintf(fileId,'%3d: ',i);
			for j=1:colSize
				fprintf(fileId,'%+6.2e ',polynomials.coef(i,j));
			end
			fprintf(fileId,'\n'); 
		end
	end
else % 1 < noOfPolynomials 
	%	fprintf(fileId,'\n***** Polynomial Systems *****\n\n'); 
    for ell=1:noOfPolynomials
        fprintf(fileId,'Polynomial %2d:\n',ell);
        fprintf(fileId,'typeCone = %+d ',polynomials{ell}.typeCone);
        fprintf(fileId,'sizeCone = %d ',polynomials{ell}.sizeCone);
        fprintf(fileId,'dimVar = %d  ',polynomials{ell}.dimVar);
        fprintf(fileId,'degree = %d  ',polynomials{ell}.degree);
        fprintf(fileId,'noTerms = %d\n',polynomials{ell}.noTerms);
        fprintf(fileId,'supportSet = \n');
        for i=1:polynomials{ell}.noTerms
            fprintf(fileId,'%3d:',i);
            for j=1:polynomials{ell}.dimVar
                fprintf(fileId,' %2d',polynomials{ell}.supports(i,j));
            end
            fprintf(fileId,'\n');
        end
        fprintf(fileId,'coefficient =\n');
        if polynomials{ell}.sizeCone == 1
            k = 0;
            for i=1:polynomials{ell}.noTerms
                k = k+1;
                if mod(k,10) == 0
                    fprintf(fileId,'\n');
                    k = 0;
                end
                fprintf(fileId,'%3d:%+6.2e',i,polynomials{ell}.coef(i,1));
            end
            fprintf(fileId,'\n');
        else
            [rowSize,colSize] = size(polynomials{ell}.coef);
            for i=1:rowSize
                fprintf(fileId,'%3d:',i);
                for j=1:colSize
                    fprintf(fileId,'%+6.2e',polynomials{ell}.coef(i,j));
                end
                fprintf(fileId,'\n');
            end
        end
    end
end
fprintf(fileId,'\n');

return

% $Header: /home/waki9/CVS_DB/SparsePOPdev/subPrograms/writeFunctions/writePolynomials.m,v 1.1.1.1 2007/01/11 11:31:51 waki9 Exp $
