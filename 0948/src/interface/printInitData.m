function printInitData(sadata, varargin)
%printInitData reports the variables (and derivatives of them) that need to 
%be initialized when solving the DAE.
%function printInitData(sadata)
%function printInitData(sadata, options)
%sadata is an object returned by daeSA.
%
%This function can have one or two optional key-value pair arguments,
%where the keys are 'varnames', and 'outfile'.
%By default, the variable names are reported as x1, x2, ...., and the
%output is on the screen
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
%'outfile', fname
%By default, the output is on the screen; if optional 'outfile',fname is
%passed, where fname is a string, the output is in file fname.
%
%Example: 
%       n = 6; G = 9.8; L = 1.0; c = 0.1;
%       sadata = daeSA(@modified2pendula,n,G,L,c);
%       cavarnames = [{'x'},{'y'},{'lam'},{'u'},{'v'},{'mu'}];
%       printInitData(sadata,'varnames',cavarnames);
%       printInitData(sadata,'varnames','z');
%
%outputs
%
% modified2pendula problem
% ----------------------------------------------------------------
% Initialization summary:
% x, x', y, y', u, v, v', v'', v'''
%
% modified2pendula problem
% ----------------------------------------------------------------
% Initialization summary:
% z1, z1', z2, z2', z4, z5, z5', z5'', z5'''
%
% See also daeSA, getInitData, printInitData
%
% Copyright 2012 Guangning Tan, Ned Nedialkov, and John Pryce

DAESAprintInitData(sadata, varargin{:});
end