function [D]=gal_diff_nonper(D,N)
% [D]=gal_diff_nonper(D,N)
% This function generates the differentiation matrix in the non-periodic
% case of Galerkin approach
% Dependencies
% dstmat_nonper.m and gal_difmatrix_nonper
D_matrix=gal_difmatrix_nonper(D,N);
C=dstmat_nonper(D,N);
D=inv(C)*D_matrix*C;
