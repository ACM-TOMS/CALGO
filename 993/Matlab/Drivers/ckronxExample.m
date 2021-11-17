
m = [40 50 60]; 
n = [30 20 10];
A = arrayfun(@(m,n) randn(m,n), m, n, 'UniformOutput', false);
B=randn(prod(n),1);
tic
C1 = ckronx(A,B);
fprintf('time using ckronx:                     %10.4f\n',toc)
tic
AA = kron(kron(A{1},A{2}),A{3});
fprintf('time to form kronecker product matrix: %10.4f\n',toc)
tic
C2 = AA*B;
fprintf('time using kronecker product matrix:   %10.4f\n',toc)


fprintf('maximum absolute difference: %10.4e\n',max(abs(C1-C2)))