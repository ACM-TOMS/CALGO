function printSolScheme(sadata, varargin)
%printSolScheme reports a solution scheme for the DAE.
%function printSolScheme(sadata)
%function printSolScheme(sadata, options)
%sadata is an object returned by daeSA.
%
%This function can have one to four optional key-value pair arguments,
%where the keys are 'varnames', 'fcnnames', 'outfile', and 'detail'.
%
%By default, the variable names are reported as x1, x2, ..., the function 
%names are reported as f1, f2, ...., and the output is on the screen.
%
%'varnames', cavarnames
%Different variable names can be given in a cell array of n strings, say
%cavarnames, and then passed as 'varnames', cavarnames. This results in the
%ith variable name being replaced by cavarnames{j}.
%
%'varnames', varname
%If each of the default names needs to be replaced by the same variable 
%name varname, followed by index number, one can pass 'varnames',varname.
%
%'fcnnames', cafcnnames
%Different function names can be given in a cell array of n strings say
%cafcnnames, and then passed as 'fcnnames', cafcnnames. This results in the
%ith function name being replaced by cafcnnames{i}.
%
%'fcnnames', fcnname
%If each of the default names needs to be replaced by the same function 
%name fcnname, followed by index number, one can pass 'fcnnames', fcnname.
%
%'outfile', fname
%By default, the output is on the screen; if optional 'outfile',fname is
%passed, where fname is a string, the output is in file fname.
%
%'detail', detaillevel
%The detaillevel can be 'compact' (default), which causes this function 
%to report the solution scheme in compact form; and 'full', which results 
%in more detailed output.
%
%  Example: 
%           n = 6; G = 9.8; L = 1.0; c = 0.1;
%           sadata = daeSA(@modified2pendula,n,G,L,c);
%           cavarnames = [{'x'},{'y'},{'lam'},{'u'},{'v'},{'mu'}];
%           cafcnnames = [{'f1'},{'g1'},{'h1'},{'f2'},{'g2'},{'h2'}];
%           printSolScheme(sadata, 'fcnnames', cafcnnames, ...
%                           'varnames',cavarnames);
%
%outputs
%
%Solution scheme for 'modified2pendula' problem
% ----------------------------------------------------------------
% Initialization summary
% x, x', y, y', u, v, v', v'', v'''
% ----------------------------------------------------------------
% k = -6: [h1] : x, y
% k = -5: [h1'] : x', y'
% k = -4: [f1, g1, h1''] : x'', y'', lam
% k = -3: [f1', g1', h1'''] : x''', y''', lam'
%         [] : v
% k = -2: [f1'', g1'', h1''''] : x'''', y'''', lam''
%         [h2] : u
%         [] : v'
% k = -1: [f1''', g1''', h1^(5)] : x^(5), y^(5), lam'''
%         [h2'] : u'
%         [] : v''
% k =  0: [f1'''', g1'''', h1^(6)] : x^(6), y^(6), lam''''
%         [h2''] : u''
%         [f2] : mu
%         [g2] : v'''
%
%See also daeSA, getInitData, getConstr, printInitData, printConstr
%
% Copyright 2012 Guangning Tan, Ned Nedialkov, and John Pryce
DAESAprintSolScheme(sadata, varargin{:});
end