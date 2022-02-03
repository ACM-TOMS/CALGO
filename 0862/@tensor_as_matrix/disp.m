function disp(t,name)
%DISP Command window display of a tensor_as_matrix.
%
%   DISP(T) displays a tensor as matrix with no name.
%
%   DISP(T,NAME) display a tensor as matrix with the given name.
%
%   Copyright 2005, Tamara Kolda and Brett Bader, Sandia National Labs
%
%   See also TENSOR_AS_MATRIX.

%Brett W. Bader and Tamara G. Kolda, Released under SAND2004-5189,
%Sandia National Laboratories, 2004.  Please address questions or
%comments to: tgkolda@sandia.gov.  Terms of use: You are free to copy,
%distribute, display, and use this work, under the following
%conditions. (1) You must give the original authors credit. (2) You may
%not use or redistribute this work for commercial purposes. (3) You may
%not alter, transform, or build upon this work. (4) For any reuse or
%distribution, you must make clear to others the license terms of this
%work. (5) Any of these conditions can be waived if you get permission
%from the authors.

if ~exist('name','var')
    name = 'ans';
end
    
namedot = [name '.'];

if strcmp(get(0,'FormatSpacing'),'compact')
    skipspaces = 1;
else
    skipspaces = 0;
end

if skipspaces ~= 1
    fprintf(1,'\n');
end

fprintf(1,'%s is a matrix corresponding to a tensor of size ',name);
printsize(t.tsize);
fprintf(1,'\n');

fprintf(1,'%s.rindices = ', name);
printvec(t.rindices);
fprintf(1,' (modes of tensor corresponding to rows)\n');
fprintf(1,'%s.cindices = ', name);
printvec(t.cindices);
fprintf(1,' (modes of tensor corresponding to columns)\n');

if skipspaces ~= 1
    fprintf(1,'\n');
end

fprintf(1,'%s',namedot);
if isempty(t.data)
    fprintf(1,'data = []\n');
else
    fprintf(1,'data = \n');
    disp(t.data);
end    

function printsize(v)
n = length(v);
for i = 1 : n - 1
    fprintf(1,'%d x ',v(i));
end
if (n > 0)
    fprintf(1,'%d', v(n));
else
    fprintf(1,'0');
end

function printvec(v)
n = length(v);
fprintf(1,'[');
for i = 1 : n - 1
    fprintf(1,'%d, ',v(i));
end
if (n > 0)
    fprintf(1,'%d', v(n));
end
fprintf(1,']');
