%% innerhs
% Discretization of the inner product in the space of weighted 
% Hilbert-Schmidt operators.
%
%% Syntax
%
%   scalProd = innerhs(A,B,seti)
%
%% Description
%
% |scalProd = innerhs(A,B,seti)| 
% computes the discretized inner product |scalProd| in the 
% space of weighted Hilbert-Schmidt operators of 
% complex matrices A, B of size seti.measNb x seti.incNb. 
% The values for the weight are stored in seti.dSMeas.
%
% * If all weights have the same value, seti.dSMeas can be a value.
% * If the weights are different, seti.dSMeas is a vector of 
%   size seti.measNb x 1. 
%   In this case A and b has to be quadratic matrices of size 
%   seti.measNb x seti.measNb.
%
% Finally, the weight |weight| is a matrix of size 
% seti.measNb x seti.measNb, in which the entries of seti.dSMeas are
% arranged on the diagonal.
%
%% Examples
%
% *Example 1: same weight value*
%
%   A = [1 4; 3 2; 0 1]; % A of size 3 x 2.
%   B = [4 3; 7 3; 1 + 2i 9]; % B of size 3 x 2.
%   seti.dSMeas = 0.1; % all weights have the same value
%   scalProd = innerhs(A,B,seti)
%
% _Result_
%
%   scalProd =
%
%       5.2000
%
% *Example 2: different weight values*
%
%   A = [1 4; 3 2]; % A of size 2 x 2.
%   B = [4 3; 1 + 2i 9]; % B of size 2 x 2.
%   seti.dSMeas = [0.1; 0.05]; % different weights
%   scalProd = innerhs(A,B,seti)
%
% _Result_
%
%   scalProd =
%
%      2.6500 - 0.3000i
%
%% Input Arguments
%
% * seti.dSMeas     : Approximation of the infinitesimal element of a closed contour 
%             with control points. It is a vector of size seti.measNb x 1.
%             If all elements are almost the same it is a value.
%             If input is only one point, it is set to 1 
%             (a closed contour does not make sense).
%             See also <expSetup.html>, <pntsGeometry.html>, and <dS2D.html>.
%
% * A               : Complex matrix of size seti.measNb x seti.incNb.
%                     (If seti.dSMeas is a vector A must be a quadratic
%                     matrix of size seti.measNb x seti.measNb).
% * B               : Complex matrix with same properties as matrix A.
%
%% Output Arguments
%
% * scalProd    : discretized inner product in the 
% space of weighted Hilbert-Schmidt operators.
%
%% More About
%
% For the space $\bf{C}^{N_s\times N_i} the discretization of the 
% inner product in space of Hilbert-Schmidt operators
% HS by weighted Frobenius product is given by
%
% $\langle A,B \rangle_\mathrm{dis} := \mathrm{trace}(B^\ast \omega^s A)$,
%
% where $N_i$ is the number of incident fields (seti.incNb),
% $N_s$ is the number of measurements (seti.measNb), and 
% $\omega^s$ are the weights (seti.dSMeas), see [1, Sec. 3.6, eq. (34)].
% 
%
%% References
%
% * [1] Florian B&uuml;rgel, Kamil S. Kazimierski, and Armin Lechleiter. A sparsity regularization and total variation based computational framework for the inverse medium problem in scattering. _Journal of Computational Physics_, 339:1-30, 2017.
%
%% See Also
% * <expSetup.html>
% * <pntsGeometry.html>
% * <dS2D.html>
% * <dS3D.html>
%
%% Code
%
function scalProd = innerhs(A,B,seti)

if length(seti.dSMeas) == 1
    scalProd = seti.dSMeas*trace(B'*A);
else
    [m,n] = size(A);
    if m ~= n
        disp('Error in innerhs.m - cannot determine correct dimension (?!)')
    end
    if m == length(seti.dSMeas)
        % Set the entries of seti.dSMeas on a matrix of size m x m
        % and store the result as sparse matrix "weight".
        weight = sparse((1:m),(1:m),seti.dSMeas);
        scalProd = trace(B'*weight*A);
    else
        disp('Error in scalProd.m - input has incorrect dimension (?!)')
    end
end
