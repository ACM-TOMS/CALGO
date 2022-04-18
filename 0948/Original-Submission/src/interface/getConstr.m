function constr = getConstr(sadata)
%getConstr returns a vector describing the constraints of
%a DAE.
%constr = getConstr(sadata) 
%sadata is an object returned by daeSA.
%
%If the DAE is structurally well-posed, getConstr returns a non-negative 
%row vector constr describing the constraints; otherwise, it returns [].
%
%Denote the i-th function in the DAE evaluation by f_i. Constraints are
%0 = f_i^(c) for c<constr(i). As a result, the sum of the elements of
%constr gives the total number of constraints.
%
%Example:
%       n = 6; G = 9.8; L = 1.0; c = 0.1;
%       sadata = daeSA(@modified2pendula,n,G,L,c);
%       constr = getConstr(sadata)
%
%produces
%       constr =
%             4     4     6    0     1     3
%
%This implies that the constraints include
%
%       4 | f1, f1', f1'', f1'''
%       4 | f2, f2', f2'', f2'''
%       6 | f3, f3', f3'', f3''', f3'''', f3^(5)
%       0 | 
%       1 | f5
%       3 | f6, f6', f6''
%
%See also daeSA, printConstr, isSWP.
%
%Copyright 2012 Guangning Tan, Ned Nedialkov, and John Pryce
constr = DAESAgetConstr(sadata);
end