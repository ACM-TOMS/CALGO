function out = iscellcomplex(M)
% ISCELLCOMPLEX(M)
% Check all entries of a cell matrix. If there is at least one complex
% number in it, then return 1. Otherwise, return 0.

out = 0;

[m,n] = size(M);
for l = 1:m
    for c = 1:n
        if sum(~isreal(M{l,c})) > 0
            out = 1;
            return
        end
    end
end

end