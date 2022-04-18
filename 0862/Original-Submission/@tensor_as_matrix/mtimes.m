function c = mtimes(a,b)
%MTIMES Multiplies two tensor_as_matrix objects.
%
%  C = MTIMES(A,B) computes the product of A and B. The result is a
%  TENSOR_AS_MATRIX object and can be transformed into a tensor.
%
%  C = MTIMES(A,B) is called for the syntax 'A * B' when A or B is a
%  TENSOR_AS_MATRIX object.  
%
%   Copyright 2005, Tamara Kolda and Brett Bader, Sandia National Labs
%
%   See also TENSOR_AS_MATRIX.

%Brett W. Bader and Tamara G. Kolda, Released under SAND2004-5189,
%Sandia National Laboratories, 2004.  Please address questions or
%comments to: tgkolda@sandia.gov.  Terms of use: You are free to copy,
%distribute, display, and use this work, under the following
%conditions. (1) You must give the original authors credit. (2) You may
%not use or redistribute this work for commercial purposes. (3) You may
%not alter, transform, or build upon this work. (4) For any reuse or
%distribution, you must make clear to others the license terms of this
%work. (5) Any of these conditions can be waived if you get permission
%from the authors.


if ~isa(a,'tensor_as_matrix')
    c = mtimes(tensor_as_matrix(a,1),b);
    return;
end

if ~isa(b,'tensor_as_matrix')
    c = mtimes(a,tensor_as_matrix(b,1));
    return;
end

if size(a,2) ~= size(b,1)
    error(['Size mismatch: Number of columns in A is not equal to' ...
	   ' the number of rows in B']);
end

c = a;
c.tsize = [a.tsize(a.rindices) b.tsize(b.cindices)];
c.rindices = 1:length(a.rindices);
c.cindices = (1:length(b.cindices)) + length(a.rindices);
c.data = a.data * b.data;

