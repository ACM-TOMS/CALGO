echo on
% ----------------------------------------------------------------------
% Example 14: Constructing a tensor by converting a
% tensor_as_matrix object.
% ----------------------------------------------------------------------
% (Note: Must run ex13 immediately before this script.)

% Convert the tensor_as_matrix object from the previous example
% into a tensor
D = tensor(B)


% They are the same
norm(C-D)

echo off
