function M = BTVMat(a)
%
% BTVMat    Returns the two 2X2 matrices introduced in [1]
%
%  M = BTVMat
%       returns M such that 
%          M{1} = [1,1;0,1]
%          M{2} = 0.78*[1,0;1,1] 
%
%  M = BTVMat(alpha)
%       returns M such that 
%          M{1} = [1,1;0,1]
%          M{2} = alpha*[1,0;1,1] 
%
% It has been shown in [1] that there are uncountably many values of alpha
% such that the family M defined above does not satisfy the finiteness 
% conjecture. 
% These matrices are also used in [2] with alpha=0.78 as an example of the 
% Balanced Complex Polytope algorithm.
%
% REFERENCES
% [1]   V.D. Blondel, J. Theys and A.A. Vladimirov, 
%         "An elementary counterexample to the finiteness conjecture",
%       SIAM J. Matrix Anal. Appl. 27(1):256-272, 2005
%
% [2]   N.Guglielmi and M.Zennaro, 
%        "Finding extremal complex polytope norms for
%         families of real matrices",
%       SIAM J. Matrix Anal. Appl. 31(2):602-620, 2009
if nargin<1
    a = 0.78;
end

M{1} = [1,1;0,1];
M{2} = a*[1,0;1,1];

end