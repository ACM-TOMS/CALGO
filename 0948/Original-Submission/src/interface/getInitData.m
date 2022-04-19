function iv = getInitData(sadata)
%getInitData returns a vector describing the variables and their 
%derivatives that need to be initialized.
%iv = getInitData(sadata)
%sadata is an object returned by daeSA.
%
%If the DAE is structurally well posed (SWP), getInitData returns a 
%non-negative row vector iv indicating the set of variables (and 
%derivatives of them) that need to be initialized when solving the DAE. 
%If the DAE is not SWP, setInitData returns [].
%
%Denote the j-th variable by x_j. iv(j) is the number of derivatives of 
%x_j that need initial values, i.e. x_j^(d) for all d<iv(j) must 
%be initialized. As a result, hte sum of the elements of iv is the total
%number of IVs needed.
%
%Example:
%       n = 6; G = 9.8; L = 1.0; c = 0.1;
%       sadata = daeSA(@modified2pendula,n,G,L,c);
%       iv = getInitData(sadata)
%
%outputs
%       iv =
%             2     2    0    1     4    0 
%
%That is, the variables that need to be initialized are
%
%       2 | x1, x1'
%       2 | x2, x2'
%       0 | 
%       1 | x4
%       4 | x5, x5', x5'', x5'''
%       0 |
%
%See also daeSA, printInitData, isSWP.
%
%Copyright 2012 Guangning Tan, Ned Nedialkov, and John Pryce

iv = DAESAgetInitData(sadata);
end