function B = diagproc(A)

%   making zero diagonal entries of A nonzero

n = size(A,2);
d = diag(A);
I = find(d);
B = A;

for j = 1:n
    if d(j) == 0
        [l,k] = min(abs(I-j));
        B(j,j) = eps*d(I(k));
    end
end
