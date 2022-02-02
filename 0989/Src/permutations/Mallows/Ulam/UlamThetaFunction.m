function[ftheta] =  UlamThetaFunction(theta,meanxjs,newton_params)
% [ftheta] =  UlamThetaFunction(theta,newton_params): Ulam Theta parameter estimation function.
%
% Input parameters:
% theta: Parameter to be evaluated
% ncount:  Vector with the number of permutations at distance d, according
%          to Ulam distance
% newton_params: Parameters used by the optimization method to find the
%                estimates
% Output parameters:
% ftheta: Value of the function
%
% References:
% [1] E. Irurozki, J. Ceberio, B. Calvo. J. A. Lozano. Sampling and learning the Mallows model under the Ulam distance. Technical report EHU-KZAA-TR;2014-04. January 2014.
% [2] J. Ceberio,  A. Mendiburu, J.A Lozano: Introducing the Mallows Model on Estimation of Distribution Algorithms. In Proceedings of International Conference on Neural Information Processing (ICONIP), 2011
%
% Created version 20/02/2015. Roberto Santana (roberto.santana@ehu.es)
%
% Last version  20/02/2015. Roberto Santana (roberto.santana@ehu.es)

n = newton_params{1};              %number of variables
ncount = newton_params{2};
distances = [1:n-1]; % 1..n-1    % Distances to the central permutation (between 1 and n-1
denom = ncount(1:n-1) ./ exp(theta.*distances); 
numer = denom .* distances;

ftheta =  sum(numer)/sum(denom) - meanxjs;


