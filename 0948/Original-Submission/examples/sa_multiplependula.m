% This is script for the multiple chained pendula example in Section 5.4 of 
% R. McKenzie, J. Pryce, N. Nedialkov, G. Tan. "DAESA User Guide"
%
% Copyright 2014 G. Tan,  N. Nedialkov,  J. Pryce 

% File: sa_multiplependula.m

p = 5;      % number of chained pendula
n = 3*p;    % size of the problem
sadata = daeSA(@multiplependula, n, p);

% display original structure
figure; showStruct(sadata);
cd Figures
print 'multiplependula_5.eps';
% display (fine) BTF
figure; showStruct(sadata,'disptype', 'fineblocks')
print 'multiplependula_5_fine.eps';
cd ..