% Sparse LDL factorization
% 
%    ldl       - LDL' factorization of a real, sparse, symmetric matrix
%    ldldemo   - demo program for LDL
%    ldlrow    - an m-file description of the algorithm used by LDL
%    ldltest   - test program for LDL
%    ldlmain2  - compiles and runs a longer test program
%
%
% LDL Version 1.1 (Apr. 22, 2005),  Copyright (c) 2005 by Timothy A. Davis.
% All Rights Reserved.
% 
% LDL License:
% 
%   Your use or distribution of LDL or any modified version of
%   LDL implies that you agree to this License.
%
%   THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
%   EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.
%
%   Permission is hereby granted to use or copy this program, provided
%   that the Copyright, this License, and the Availability of the original
%   version is retained on all copies.  User documentation of any code that
%   uses LDL or any modified version of LDL code must cite the
%   Copyright, this License, the Availability note, and "Used by permission."
%   Permission to modify the code and to distribute modified code is granted,
%   provided the Copyright, this License, and the Availability note are
%   retained, and a notice that the code was modified is included.  This
%   software was developed with support from the National Science Foundation,
%   and is provided to you free of charge.
%
% Availability:
%
%   http://www.cise.ufl.edu/research/sparse/ldl
%
% Acknowledgements:
%
%   This work was supported by the National Science Foundation, under
%   grant CCR-0203270.
%
%   Portions of this work were done while on sabbatical at Stanford University
%   and Lawrence Berkeley National Laboratory (with funding from the SciDAC
%   program).  I would like to thank Gene Golub, Esmond Ng, and Horst Simon
%   for making this sabbatical possible.
