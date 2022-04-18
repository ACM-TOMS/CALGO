function printConstr(sadata, varargin)
%printConstr reports the functions (and derivatives of them) that form the 
%constraints.
%function printConstr(sadata)
%function printConstr(sadata, options)
%sadata is an object returned by daeSA.
%
%This function can have one or two optional key-value pair arguments,
%where the keys are 'fcnnames' and 'outfile'.
%By default, the function names are reported as f1, f2, ...., and the 
%output is on the screen. 
%
%'fcnnames', cafcnnames
%Different function names can be given in a cell array of n strings, say
%cafcnnames, and then passed as 'fcnnames', cafcnnames. This results in the
%ith function name being replaced by cafcnnames{i}.
%
%'fcnnames', fcnname
%If each of the default names needs to be replaced by the same function 
%name fcnname, followed by index number, one can pass 'fcnnames', fcnname.
%
%'outfile', fname
%By default, the output is on the screen. If optional 'outfile',fname
%is passed, where fname is a string, the output is stored in file fname.
%
%Example: 
%       n = 6; G = 9.8; L = 1.0; c = 0.1;
%       sadata = daeSA(@modified2pendula,n,G,L,c);
%       cafcnnames = [{'f1'},{'g1'},{'h1'},{'f2'},{'g2'},{'h2'}];
%       printConstr(sadata,'fcnnames',cafcnnames);
%       printConstr(sadata,'fcnnames','g');
%
%outputs
%
% modified2pendula problem
% ---------------------------------------------------------------------------
% Constraints:
% f1, f1', f1'', f1''', g1, g1', g1'', g1''', h1, h1', h1'', h1''', h1'''', 
% h1^(5), g2, h2, h2', h2''
%
% modified2pendula problem
% ---------------------------------------------------------------------------
% Constraints:
% g1, g1', g1'', g1''', g2, g2', g2'', g2''', g3, g3', g3'', g3''', g3'''', 
% g3^(5), g5, g6, g6', g6''
%
% See also daeSA, getConstr, printSolScheme
%
% Copyright 2012 Guangning Tan, Ned Nedialkov, and John Pryce

DAESAprintConstr(sadata, varargin{:});
end