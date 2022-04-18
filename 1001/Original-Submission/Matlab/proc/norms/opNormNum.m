%% opNormNum
%
% Numerical approximation of the operator norm L = ||K||.
%
%% Syntax
%
%   L = opNormNum(xnRVD,JA,JB,KcompNorm,seti,dispDepth)
%
%% Description
%
% |L = opNormNum(xnRVD,JA,JB,KcompNorm,seti)| computes the operator norm
% $L = \|K\|$ approximately and multiplied by safety factor 2.
%
% This function is called in <pda.html>
%
%% Input Arguments
%
% * xnRVD   :   solution in <pda.html> 
%               (real vector of size |2*seti.nInv^seti.dim x 1|)
%               (R: real, V: vector: D: scaled down ROI).
% * JA, JB  :   Jacobi matrices, see <mimo.html>
% * KcompNorm   :   numerator in operator norm definition; 
%                   is prepared in subfunction |KcomponentsStruct| in <pda.html>.
% * seti        :   structural array
% * seti.dVinv  :   is needed in <innerinv.html>, called in <norminv2.html>
% * dispDepth   :   controlls the depth of displayed messages.
%
%% Output Arguments
%
% * L   :   approximately computed operator norm multiplied by safety factor 2
%
%% More About
%
% *Theory: operator norm L*
%
% The operator norm L is essentially
%
% $$ L = \|K\| = \max_{x \neq 0} 1/\|x\|_2 \ \sqrt{\sum_i \|K_i(x)\|_2^2} $$
%
% with the components $K_i$ of $K: x \to Y$.
%
% Note that the norms must be induced by inner products, see [1].
%
% The components of $K$ was already stored in 
%
% $\texttt{KcompNorm} =  \sqrt{\sum_i \|K_i(x)\|_2^2}$.
%
% Of course in nominator and denominator can appear different norms: the
% notation $\|\cdot\|_2$ was choosen to emphasize the requirement of a norm
% induced by a inner product.
%
% In the case of |seti.invNo = 6|, which is the default case in the 
% published version, see <setInvType.html>, 
%
% $L = \|K\| =
%  \max_{x \neq 0} \sqrt{\|K_{\mathrm{dis}}(x)\|_2^2 + \|K_{\mathrm{tv}}(x)\|_2^2} / \|x\|_2$
%
% with norms <normws2.html>, <normTVinv2.html>, and <norminv2.html> for
% $\|\cdot\|_2$ (in this order).
%
% See also [2, Sec. 4.8].
%
% *Numerics: operator norm L*
%
% * To compute an approximation for L we test the operator with two kinds
% of sample vectors:
% * Sample vectors 1: random vectors
% * Sample vectors 2: vectors of the form (1, . . . , 1, 0, . . . , 0).
% 
% * Note that the approximately computed operator norm is too small in most cases...
% * Then stepsize in pda algorithm is too big (that is a problem).
% * To avoid too big stepsize in pda algorithm we use safety factor 2
%
%% References
%
% * [1] Antonin Chambolle and Thomas Pock. A first-order primal-dual algorithm for convex problems with applications to imaging. _Journal of Mathematical Imaging and Vision_, 40(1):120-145, 2011.
% * [2] Florian B&uuml;rgel, Kamil S. Kazimierski, and Armin Lechleiter. A sparsity regularization and total variation based computational framework for the inverse medium problem in scattering. _Journal of Computational Physics_, 339:1-30, 2017.
%
%% See Also
%
% * <pda.html>
% * <mimo.html>
% * <norminv2.html>
% * <innerinv.html>
%
%% Code

function L = opNormNum(xnRVD,JA,JB,KcompNorm,seti,dispDepth)

%%
% *Random vectors*
iN = 100;
opnorm1 = zeros(iN,1);
for i = 1:iN
    h = rand(size(xnRVD))+1i*rand(size(xnRVD));
    opnorm1(i) = KcompNorm(h,JA,JB)/norminv2(h,seti);
end
L1 = max(opnorm1);

%%
% *Vectors with one and zero...*

% iN = length(x0); % too expensive
iN = 1;
opnorm2 = zeros(iN,1);
h = zeros(size(xnRVD));
for i = 1:iN % x = [1, 0, ..., 0], x = [1, ..., 1, 0, ..., 0] until x = [1, ..., 1]
    h(i) = 1;
    opnorm2(i) = KcompNorm(h,JA,JB)/norminv2(h,seti);
end
L2 = max(opnorm2);

%%
% *Final operator norm*
opnorm = max(L1,L2);
L = 2*max(opnorm); % multiplication with safety factor 2
if dispDepth >= 3
    disp('      Operator norm L computed approximately and multiplied by safety factor 2.')
end
end
