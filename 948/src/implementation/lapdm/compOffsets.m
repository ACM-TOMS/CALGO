function [c, d, ccl, dcl, HVT] = compOffsets(A, hvt, r)

n_blocks = length(r)-1;
n = size(A, 1);
c = zeros(n, 1); d = zeros(1, n); % global offsets
ccl = zeros(n, 1); dcl = zeros(1, n); % global offsets
HVT = zeros(1, n);

% Compute offsets block by block
for block = 1:n_blocks
    crows = 1:r(block)-1;
    ccols = r(block):r(block+1)-1;
    
    ncols = r(block+1) - r(block); % Size of block
    
    D = A(ccols, ccols);
    if block>1
        C = A(crows, ccols) + repmat(c(crows), [1 ncols]);
        d_init1 = max([C; D], [], 1);
        d_init2 = max(D, [], 1);
        % Find global offsets
        [c(ccols), d(ccols)] = findcd(D, hvt(ccols), d_init1);
        % Find coarse local offsets
        [ccl(ccols), dcl(ccols)] = findcd(D, hvt(ccols), d_init2);
    else % first block
        d_init1 = max(D, [], 1);
        [c(ccols), d(ccols)] = findcd(D, hvt(ccols), d_init1);
        ccl(ccols) = c(ccols); dcl(ccols) = d(ccols);
    end
    
    % Permute HVT
    HVT(ccols) = hvt(ccols) + r(block) - r(1);
end
end