'testProjOntoPSDcone'


%% not enough input
X = []; PSDcone = []; rankUbd = [];
exception = [];
try
    [Xp, rank] = projOntoPSDcone;
catch exception
    display(['Catched:' exception.message]);
end
if isempty(exception); error('Uncatched: not enough input'); end


exception = [];
try
    [Xp, rank] = projOntoPSDcone(X);
catch exception
    display(['Catched:' exception.message]);
end
if isempty(exception); error('Uncatched: not enough input'); end

%% not enough output
exception = [];
try
    projOntoPSDcone(X, PSDcone);
catch exception
    display(['Catched:' exception.message]);
end
if isempty(exception); error('Uncatched: not enough output'); end

%% inconsistent inputs
X = zeros(3^2 + 5^2, 1);
PSDcone = [3 4];
try
    Xp = projOntoPSDcone(X, PSDcone);
catch exception
    display(['Catched:' exception.message]);
end
if isempty(exception); error('Uncatched: inconsistent inputs'); end

PSDcone = [3 5];
rankUbd = [2 2 2];
try
    Xp = projOntoPSDcone(X, PSDcone, rankUbd);
catch exception
    display(['Catched:' exception.message]);
end
if isempty(exception); error('Uncatched: inconsistent inputs'); end

%% empty output
PSDcone = [10000 10000];
X = zeros(0, 0); 
Xpa = projOntoPSDcone(X, PSDcone);
[XpA, rankA] = projOntoPSDcone(X, PSDcone);

X = zeros(10, 0); 
Xpb = projOntoPSDcone(X, PSDcone);
[XpB, rankB] = projOntoPSDcone(X, PSDcone, PSDcone);

X = zeros(0, 10); 
Xpc = projOntoPSDcone(X, PSDcone);
[XpC, rankC] = projOntoPSDcone(X, PSDcone, PSDcone);
if isempty(Xpa) && isempty(XpA) && isempty(rankA) && ...
        isempty(Xpb) && isempty(XpB) && isempty(rankB) && ...
        isempty(Xpc) && isempty(XpC) && isempty(rankC) 
    disp('empty output: OK');
else
    error('Error: non empty output'); 
end

%% already PSD
X1 = eye(3); X2 = eye(5);
X = [X1(:); X2(:)];
PSDcone = [3, 5];
Xp = projOntoPSDcone(X, PSDcone);
if norm(Xp - X) < eps
    disp('Positive Definite: OK');
else
    error('Positive Definite: Error');
end

X(end) = 0;
[Xp, rank] = projOntoPSDcone(X, PSDcone);
if norm(Xp - X) < eps && rank(1) == 3 && rank(2) == 4
    disp('Positive Semi-Definite: OK');
else
    error('Positive Semi-Definite: Error');
end

rng(1)
n = 10;
x = randn(n, 1) * 2;
X = x * x';
PSDcone = n;
[Xp, rank] = projOntoPSDcone(X, PSDcone);
if norm(Xp(:) - X(:)) < eps * norm(X(:)) * 10 && rank == 1
    disp('Rank-1 PSD: OK');
else
    error('Rank-1 PSD: error');
end

%%
n = 10;
X = (rand(n, 1) - 0.5) * 20;
PSDcone = ones(n, 1);
[Xp, rank] = projOntoPSDcone(X, PSDcone);
if norm(Xp - X.* (X>0)) < eps * norm(X(:)) && isequal(rank, X>0)
    disp('Scalar blocks: OK');
else
    error('Scalar blocks: error');
end

