function X = genPerms(n, deg, homogeneous)

% GENPERMS Generation of all non-negative integer n-tuples of given (max) sum
%
%    X = genPerms(n, deg) returns a matrix containing all non-negative integer
%      n-tuples such that the sum of their components is <= deg. The n-tuples
%      will be lexicographically ordered.
%
%    X = genPerms(n, deg, 1) returns a matrix containing all non-negative
%      integer (n+1)-tuples such that the sum of their components is == deg.
%      The (n+1)-tuples will be lexicographically ordered.
%

X = generate(zeros(nck(n + deg, n), n), 1, deg, n, 1, zeros(1, n));

if (nargin == 3) && (homogeneous == 1),
    X = [X, deg-sum(X, 2)];
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X, i] = generate(X, i, deg_left, nb_var, cur, prefix)

if cur > nb_var,
    X(i, :) = prefix;
    i = i + 1;
else
    for k=0:deg_left,
        prefix(cur) = k;
        [X, i] = generate(X, i, deg_left - k, nb_var, cur + 1, prefix);
    end
end
end