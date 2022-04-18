% Tensor Toolbox.
% Efficient manipulation of multidimensional arrays for Matlab.
% Version 1.0  13-April-2006
% by Tamara G. Kolda and Brett W. Bader
% Copyright 2005, Sandia National Labs
% 
% @TENSOR
%   and          - Logical AND for tensors.
%   disp         - Command window display of a tensor.
%   display      - Command window display of a tensor.
%   double       - Convert a tensor to double array.
%   end          - Last index of indexing expression for tensor.
%   eq           - Equal for tensors.
%   ge           - Greater than or equal for tensors.
%   gt           - Greater than for tensors.
%   issamesize   - Compare sizes for tensors.
%   ldivide      - Left array divide for tensor.
%   le           - Less than or equal for tensor.
%   lt           - Less than for tensor.
%   minus        - Binary subtraction for tensors.
%   mtimes       - Scalar or inner product for tensors.
%   multiarrayop - Generic multidimensional array functions for tensors. 
%   ndims        - Return the number of dimensions of a tensor.
%   ne           - Not equal for tensors.
%   norm         - Frobenius norm of a tensor.
%   not          - Logical NOT for tensors.
%   or           - Logical OR for tensors.
%   order        - Return the order of a tensor.
%   permute      - Permute tensor dimensions.
%   plus         - Binary addition for tensors. 
%   power        - Elementwise power operator for a tensor.
%   rdivide      - Right array divide for tensors.
%   shiftdim     - Shift dimensions for a tensor.
%   size         - Size of tensor.
%   squeeze      - Remove singleton dimensions from a tensor.
%   subsasgn     - Subscripted assignment for a tensor.
%   subsref      - Subscripted reference for tensors.
%   tensor       - Tensor object constructor.
%   times        - Element-wise multiplication for tensors.
%   ttm          - Tensor times matrix.
%   ttt          - Tensor mulitplication (tensor times tensor).
%   ttv          - Tensor times vector.
%   uminus       - Unary minus for tensors.
%   uplus        - Unary plus for tensors.
%   xor          - Logical EXCLUSIVE OR.
%
% @TENSOR_AS_MATRIX
%   ctranspose       - Complex conjugate transpose for tensor_as_matrix.
%   disp             - Command window display of a tensor_as_matrix.
%   display          - Command window display of a tensor_as_matrix.
%   double           - Convert tensor_as_matrix to double array.
%   mtimes           - Multiplies two tensor_as_matrix objects.
%   size             - Size of tensor_as_matrix.
%   subsasgn         - Subscripted assignment for tensor_as_matrix.  
%   subsref          - Subscripted reference for tensor_as_matrix.
%   tensor_as_matrix - Constructor for matrix representation of a tensor
%   tsize            - Tensor size of tensor_as_matrix.
%
% @CP_TENSOR
%   cp_tensor  - Tensor stored in CANDECOMP/PARAFAC form.
%   disp       - Command window display for a cp_tensor.
%   display    - Command window display for a cp_tensor.
%   full       - Convert a cp_tensor to a (dense) tensor.
%   issamesize - Size comparison for cp_tensor.
%   minus      - Binary subtraction for cp_tensor.  
%   mtimes     - Implement A*B (scalar multiply) for cp_tensor.
%   ndims      - Number of dimensions for a cp_tensor.
%   order      - Number of dimensions for a cp_tensor.
%   permute    - Permute dimensions of a cp_tensor.
%   plus       - Binary addition for cp_tensor.
%   size       - Size of cp_tensor.
%   subsasgn   - Subscripted assignement for cp_tensor.
%   subsref    - Subscripted reference for a cp_tensor.
%   ttv        - Tensor times vector for cp_tensor.
%   uminus     - Unary minus for cp_tensor. 
%   uplus      - Unary plus for a cp_tensor. 
%   arrange    - Arranges the rank-1 terms of a CP tensor.
%   norm       - Frobenius norm of a CP tensor.
%
% @TUCKER_TENSOR
%   disp          - Command window display of a tucker_tensor
%   display       - Command window display of a tucker_tensor.
%   full          - Convert a tucker_tensor to a (dense) tensor.
%   issamesize    - Comparison of sizes for a tucker_tensor and another object.
%   mtimes        - Implement scalar multiplication for a tucker_tensor.
%   ndims         - Return the number of dimensions for a tucker_tensor.
%   order         - Return the number of dimensions for a tucker_tensor.
%   permute       - Permute dimensions for a tucker_tensor.
%   size          - Size of a tucker_tensor.
%   subsasgn      - Subscripted reference for a tucker_tensor.
%   subsref       - Subscripted reference for a tucker_tensor.
%   tucker_tensor - Tensor stored in Tucker form.
%   uminus        - Unary minus for tucker_tensor.
%   uplus         - Unary plus for tucker_tensor.
%
% EXTRAS
%   khatrirao    - Khatri-Rao product
%   revkhatrirao - Reverse Khatri-Rao product
%
% EXAMPLES
%   ex1,..,ex23  - Demos
%   setup        - Set-up for demos
%   runall       - Run all demos
%   hopm         - Higher-order power method
%   hooi         - Higher-order orthogonal iteration

