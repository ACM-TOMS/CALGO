function disp(t, name)
%DISP Command window display for a cp_tensor.
%
%   DISP(T) displays a CP tensor with no name.
%
%   DISP(T,NAME) display a CP tensor with the given name.
%
%   Copyright 2005, Tamara Kolda and Brett Bader, Sandia National Labs
%
%   See also DISP, CP_TENSOR/DISPLAY, CP_TENSOR

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

% Check to see whether the tensor should be printed in compact mode
compact = strcmp(get(0,'FormatSpacing'),'compact');

if ~compact
    fprintf(1,'\n');
end

fprintf(1,'%s is a CP tensor of size ', name);
printsize(size(t));

if ~compact
    fprintf(1,'\n');
end


disp(' ');
disp([name, '.lambda = ']);
disp(t.lambda);

for j = 1 : order(t)
    disp([name, '.U{', int2str(j), '} = ']);
    disp(t.u{j});
end


%----------------------------------------
function printsize(sz)

for i = 1 : length(sz) - 1
    fprintf(1,'%d x ',sz(i));
end
fprintf(1,'%d', sz(length(sz)));

