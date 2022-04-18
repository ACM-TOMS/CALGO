function [x, k] = deperiod(X)
%
% DEPERIOD  finds a periodic sequence of indexes
%
%   [Y, K] = DEPERIOD(X)
%      Where X is a vector.
%      If X contains a periodic sequence, it is returned in Y 
%      and its length is returned in K.
%      If X is not periodic, Y = X and K = length(X).
%      
% ---------------------------
% Example :
%    x = [1 2 1 2 1 2];
%    [y, k] = deperiod(x)
%
% returns
%   y = [1 2]
%   k = 2
% ---------------------------
%

n = length(X);

for k = 1:n/2,
    if (mod(n, k) == 0),
        ok = 1;
        for t = 1:n-k,
            if X(t) ~= X(t+k),
                ok = 0;
                break;
            end
        end
        if ok,
            x = X(1:k);
            return;
        end
    end
end

x = X;
k = n;

end
