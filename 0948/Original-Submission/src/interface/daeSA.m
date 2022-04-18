function sadata = daeSA(daefcn, n, varargin)
%daeSA performs structural analysis of a DAE.
%sadata = daeSA(daefcn, n, varargin)
%daefcn is a function for evaluating the DAE, n is the problem size, and 
%varargin is an optional list of one or more parameters that are passed to
%daeSA.
%daeSA returns an object that contains data from the structural analysis.
%This object is an input parameter to the remaining functions of DAESA; 
%see the list below.
%
%Example:
%       n = 6; G = 9.8; L = 1.0; c = 0.1;
%       sadata = daeSA(@modified2pendula,n,G,L,c);
%
%See also
%daeSA, getConstr, getDAEfhandle, getDOF, getIndex, getInitData, 
%getMissing, getOffsets, getBTF, getQLdata, getSigma, getSize,
%isSWP, printConstr, printInitData, printSolScheme, showStruct
%
% Copyright 2012 Guangning Tan, Ned Nedialkov, and John Pryce

sadata = DAESAmain(daefcn, n, varargin{:});
end
