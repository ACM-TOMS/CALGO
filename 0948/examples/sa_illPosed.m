% This script shows the structure of the two ill-posed problems
% discussed in Section 5.2 of R. McKenzie, J. Pryce, N. Nedialkov, G. Tan. 
% "DAESA User Guide"
%
% Copyright 2014 G. Tan, N. Nedialkov, J. Pryce

% File: sa_illPosed.m

n = 6; G = 9.8; L = 1.0; c = 0.1;

% missing equation
figure(11);
sadata1 = daeSA(@illPosed1,n,G,L,c);
showStruct(sadata1);
cd Figures
print 'illPosed1.eps';

% missing variable
figure(12);
sadata2 = daeSA(@illPosed2,n,G,L,c);
showStruct(sadata2);
print 'illPosed2.eps';

format compact; 
% has under- and over-determined parts
figure(13);
sadata3 = daeSA(@illPosed3,6);
[pe,pv,cdb,fdb] = getBTF(sadata3)
showStruct(sadata3,'disptype','blocks');
print 'illPosed3.eps';
cd ..