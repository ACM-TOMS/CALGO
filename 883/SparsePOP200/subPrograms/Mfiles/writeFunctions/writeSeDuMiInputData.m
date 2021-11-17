function writeSeDuMiInputData(fileId,printLevel,A,b,c,K)
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

fprintf(fileId,'# SeDuMi format SDP SDP data\n'); 

[m,n] = size(A); 
ell = length(K.s);
if isfield(K,'f') 
	fprintf(fileId,'K.f = %3d\n',K.f);
end
if isfield(K,'l') 
	fprintf(fileId,'K.l = %3d\n',K.l);
end
fprintf(fileId,'K.s = ');
for k=1:ell
	fprintf(fileId,'%3d',K.s(k));
end
fprintf(fileId,'\n'); 
if printLevel >= 2 
	fprintf(fileId,'Objective coefficient vector b to be maximized = \n'); 
	jj = 0;
	for j=1:m
		fprintf(fileId,'%2d:%+6.2f ',j,b(j,1));
		jj=jj+1; 
		if mod(jj,10) == 0
			fprintf(fileId,'\n');
			jj=0;
		end
	end
	fprintf(fileId,'\n');
else 
	fprintf(fileId,'Dimension of the objective coefficient vector b = %d\n',length(b'));
end
if printLevel >= 2
	fprintf(fileId,'Right hand constant vector c = \n'); 
	jj = 0; 
	for i=1:n
		fprintf(fileId,'%2d:%+6.2f ',i,c(i,1));
		jj=jj+1; 
		if mod(jj,10) == 0
			fprintf(fileId,'\n');
			jj=0;
		end
	end
	fprintf(fileId,'\n');
else
	fprintf(fileId,'Dimension of the right hand constant vector c = %d\n',length(c')); 
end
if printLevel >= 3 
	fprintf(fileId,'Coefficient matrix A = \n   '); 
	for j=1:n
		fprintf(fileId,'%2d',j);
	end
	fprintf(fileId,'\n');
	for i=1:m
		fprintf(fileId,'%2d:',i); 
		for j=1:n
			fprintf(fileId,'%+2d',A(i,j));
		end
		fprintf(fileId,'\n');
	end
else
	fprintf(fileId,'The size of the coefficient matrix A = (%d,%d)\n',m,n); 
end
fprintf(fileId,'\n'); 
return



% $Header: /home/waki9/CVS_DB/SparsePOPdev/subPrograms/writeFunctions/writeSeDuMiInputData.m,v 1.1.1.1 2007/01/11 11:31:51 waki9 Exp $
