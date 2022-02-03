function poly = toSPOPformatPoly(poly, typeCone)
%%TOSparsePOPformat Converts a poly described in the simplified SparsePOP format 
%%to the poly described in the SparsePOP format.
%%
%% It computes the fields 'noTerms', 'dimVar', and 'degree' using 'supports'.
%%
%% Usage:
%%   poly = toSparsePOPformat(poly); % set typeCone to -1.
%%   poly = toSparsePOPformat(poly, typeCone);
%%
%% Simplified SparsePOP format
%% poly.supports = a set of supports of f (x),
%%                 a poly.noTerms Å~ poly.dimVar matrix.
%% poly.coef     = coefficients,
%% 
%% SparsePOP format
%% poly.typeCone = 1 if f(x) \in R[x] is used as an objective function,
%%               = 1 if f(x) \in R[x] is used as an inequality constraint f(x) ? 0, 
%%               =-1 if f(x) \in R[x] is used as an equality constraint f(x) = 0. 
%% poly.sizeCone = 1 if poly.coef is a column vector with dimensiion poly.noTerms,
%%                 m if poly.coef is a poly.noTerms x m matrix. 
%% poly.degree   = the degree of f (x).  
%% poly.dimVar   = the dimension of the variable vector x.
%% poly.noTerms  = the number of terms of f (x).
%% poly.supports = a set of supports of f (x),
%%                 a poly.noTerms Å~ poly.dimVar matrix.
%% poly.coef     = coefficients,
%%                 a column vector of poly.noTerms dimension.
    if nargin == 2
        poly.typeCone = typeCone;
    end
    poly.sizeCone = 1; % size(poly.coef,2);
    poly.noTerms = size(poly.supports, 1);
    poly.dimVar = size(poly.supports, 2);
    poly.degree = full(max(sum(poly.supports, 2)));
end