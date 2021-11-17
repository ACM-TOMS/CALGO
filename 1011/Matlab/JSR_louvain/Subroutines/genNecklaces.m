function X = genNecklaces(n, k)

% GENNECKLACES Generation of all n-bead necklaces with k colors
%
%    X = genNecklaces(n) returns a matrix containing all n-bead necklaces
%      with 2 colors, i.e., all equivalence classes of binary n-tuples under
%      cyclic permutation (rotation), each equivalence class being represented
%      by the lexicographically first tuple in the class.
%      If n is a vector, the computation will be done on each value given in n,
%      and the returned value will be a cell array of matrices.
%
%    X = genNecklaces(n, k) does the same with k-array n-tuples.
%

% Default argument
if (nargin == 1),
    k = 2;
end
w = length(n);
X = cell(w, 1);

% Brute-force generation
for cur = 1:w,
    N = n(cur);
    scaler = k.^(N-1:-1:0)';
    % Generate all k-ary N-tuples
    A = zeros(k^N, N);
    for idx = 1:N,
        w = k^(N-idx);
        for i = 0:(k^idx)-1,
            A(i*w+1:(i+1)*w, idx) = mod(i, k);
        end
    end
    % Remove all non-necklaces
    ok = true(k^N, 1);
    shifter = zeros(N-1, N);
    for i = 1:N-1,
        shifter(i, :) = [(i+1):N, 1:i];
    end
    for scan = 1:k^N,
        if ok(scan),
            base = A(scan, :);
            ok(1+base(shifter)*scaler) = 0;
            ok(scan) = 1;
        end
    end
    X{cur} = A(ok, :);

end

if (w == 1),
    X = X{1};
end
end
