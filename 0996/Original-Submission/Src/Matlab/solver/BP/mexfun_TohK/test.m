rng(1)
numblk = 3;
PSDcone = [5 4 3];
Acell = cell(numblk, 1);
Vcell = cell(numblk, 1);
Dcell = cell(numblk, 1);
for kk = 1:numblk
    Ak = randn(PSDcone(kk));
    Acell{kk} = Ak + Ak';
    [Vcell{kk}, Dcell{kk}] = eig(Acell{kk});
end
Vvec2 = cellfun(@(x) x(:), Vcell, 'UniformOutput', false);
Vvec2 = cell2mat(Vvec2);
d2 = cellfun(@(x) diag(x), Dcell, 'UniformOutput', false);
d2 = cell2mat(d2);

Avec = cellfun(@(x) x(:), Acell, 'UniformOutput', false);
Avec = cell2mat(Avec);
[Vvec, d] = blkEigMex(Avec, PSDcone);

if norm(Vvec - Vvec2) > eps
    fprintf('error1');
end


%%
PSDcone = [5 4 0 3]; 
[Vvec, d] = blkEigMex(Avec, PSDcone);

if norm(Vvec - Vvec2) > eps
    fprintf('error2');
end

%%
Avec = randn(10, 0);
PSDcone = 0;
[Vvec, d] = blkEigMex(Avec, PSDcone);

PSDcone = [];
[Vvec, d] = blkEigMex(Avec, PSDcone);

%%
Avec = randn(0, 10);
PSDcone = 0;
[Vvec, d] = blkEigMex(Avec, PSDcone);

PSDcone = [];
[Vvec, d] = blkEigMex(Avec, PSDcone);

%%
Avec = [];
PSDcone = 0;
[Vvec, d] = blkEigMex(Avec, PSDcone);

%%
K1.s = PSDcone;
Xp = blkprojSDPBP(Avec, K1);

Xp2 = myblkprojSDPBP(Avec, K1);
norm(Xp-Xp2)