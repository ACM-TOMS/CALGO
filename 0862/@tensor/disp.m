function disp(t,name)
%DISP Command window display of a tensor.
%
%   DISP(T) displays a tensor with no name.
%
%   DISP(T,NAME) displays a tensor with the given name.
%
%   Copyright 2005, Tamara Kolda and Brett Bader, Sandia National Labs
%
%   See also TENSOR, TENSOR/DISPLAY.

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
name = [name ' is a tensor'];

% Check to see whether the tensor should be printed in compact mode
compact = strcmp(get(0,'FormatSpacing'),'compact');

if ~compact
    fprintf(1,'\n');
end

fprintf(1,'%s of ',name);
printsize(t.size);
fprintf(1,'\n');

if ~compact
    fprintf(1,'\n');
end

fprintf(1,'%s',namedot);
if isempty(t.data)
    fprintf(1,'data = []\n');
else
    fprintf(1,'data = \n');
    disp(t.data);
end    

function printsize(sz)

if isempty(sz)
    fprintf(1,'order 0 (i.e., a scalar)');
    return;
end

fprintf(1,'size ');
for i = 1 : length(sz) - 1
    fprintf(1,'%d x ',sz(i));
end
fprintf(1,'%d', sz(length(sz)));

