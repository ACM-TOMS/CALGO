% This script shows the structure of the chemical Akzo Nobel problem
% discussed in Section 5.3 of R. McKenzie, J. Pryce, N. Nedialkov, G. Tan. 
% "DAESA User Guide"
%
% Copyright 2014 G. Tan, N. Nedialkov, J. Pryce

% File: sa_chemakzo.m

figure(11);
sadata = daeSA(@akzonobel,6);
showStruct(sadata,'disptype','fineblocks');
cd Figures
print 'chemakzo_finebloks.eps';

cd ../Reports
if (exist('chemakzo.txt', 'file'))
    delete('chemakzo.txt');
end;
printSolScheme(sadata,'outfile','chemakzo.txt','varnames','y');
cd ..