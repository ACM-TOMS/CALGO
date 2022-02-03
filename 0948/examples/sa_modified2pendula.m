% This is script for the modified coupled pendula example in Section 5.1 of 
% R. McKenzie, J. Pryce, N. Nedialkov, G. Tan. "DAESA User Guide"
%
% Copyright 2014 G. Tan,  N. Nedialkov,  J. Pryce 

% File: sa_modified2pendula.m

% do structural analysis
%{startsadata}
n = 6; G = 9.8; L = 1.0; c = 0.1;
sadata = daeSA(@modified2pendula,n,G,L,c);
%{endsadata}
% display original structure
figure(11);
%{startshowmod2pend}
showStruct(sadata);
%{endshowmod2pend}
cd Figures
print 'modified2pendula.eps';
% display coarse blocks
figure(12);
%{startcoarsemod2pend}
showStruct(sadata,'disptype','blocks');
%{endcoarsemod2pend}
print 'modified2pendulaDMcoarse.eps';

% display fine blocks
figure(13);
%{startfinemod2pend}
showStruct(sadata,'disptype','fineblocks');
%{endfinemod2pend}
print 'modified2pendulaDMfine.eps';

% report index and DOF
%{startindexandDOF}
index = getIndex(sadata);
DOF = getDOF(sadata);
%{endindexandDOF}

% report variables that need to be initialized 
cd ../Reports
if (exist('initdata.txt', 'file'))
    delete('initdata.txt');
end;
%{startprintInitData}
vars = {'x','y','lam','u','v','mu'};
printInitData(sadata,'varnames',vars,...
                'outfile','initdata.txt');
%{endprintInitData}
			
% report constraints
if (exist('constr.txt', 'file'))
    delete('constr.txt');
end;
printConstr(sadata,'outfile','constr.txt');

% report solution scheme
if (exist('solscheme.txt', 'file'))
    delete('solscheme.txt');
end;
%{startprintSolScheme}
printSolScheme(sadata,'varnames',vars,...
                      'outfile','solscheme.txt','detail','full');
%{endprintSolScheme}
        
% report solution scheme in compact form
if (exist('solschemecompact.txt', 'file'))
    delete('solschemecompact.txt');
end;
%{startprintSolSchemeCompact}
printSolScheme(sadata,'varnames',vars,...
               'outfile','solschemecompact.txt','detail','compact');
%{endprintSolSchemeCompact}
cd ..