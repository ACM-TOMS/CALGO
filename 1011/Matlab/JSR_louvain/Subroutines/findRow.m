function [idx, insert] = findRow(M, row)

% FINDROW Finds a given row in an ordered matrix.
%
%    idx = findRow(M, r) returns an index idx such that M(idx, :) == r, if it
%      is possible, otherwise idx == 0. The rows of the matrix M are supposed
%      to be lexicographically ordered.
%
%    [idx, insert] = findRow(M, r) moreover returns the index of the row where
%      r should be inserted if it is not in M. Otherwise, insert == 0.
%

if (isempty(M)),
    idx = 0;
    insert = 1;
else
    m = size(M, 2);
    low = 1;
    high = size(M, 1) + 1;
    
    % bisection
    while (low ~= high),
        mid = floor((low + high) / 2);
        
        for i=1:m+1,
            if (i == m + 1),
                idx = mid;
                insert = 0;
                return;
            end
            if row(i) < M(mid, i),
                high = mid;
                break;
            elseif row(i) > M(mid, i),
                low = mid + 1;
                break;
            end
        end
    end
    
    % Line not found
    idx = 0;
    insert = low;
end
