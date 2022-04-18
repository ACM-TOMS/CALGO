% EXPERIMENT2 - Produce results given in Table 2.
%
% DESCRIPTION:
%   Compare phipm and phip on various large sparse matrices.

% ===========================================================
NumExp = 100;
p = 3;
% ===========================================================
% 6.2 - A nonsymmetric example
fprintf('*** 6.2 - A nonsymmetric example orani678 ***\n');
P62 = load('orani678')
T0 = 0;
TF = 10;
N = 2529;

% EXACT
fprintf('------------------------------\n');
fprintf('*** Exact ***\n');
K = spdiags(ones(p, 1), 1, p, p);
A = [P62.Problem.A, ones(N, p); zeros(p, N), K];
ep = zeros(p,1); ep(end) = 1;
b = [ones(2529, 1); ep];
[wexact2, stats] = phipm(TF, P62.Problem.A, ones(N, p+1), eps, 0);
[wexact, tvals, mv] = expmv_tspan(A, b, T0, TF, 100, 'double');
wexact = wexact(1:N,end);
norm(wexact2-wexact)

% PHIP
fprintf('------------------------------\n');
fprintf('*** PHIP ***\n');
for i = 1:NumExp
  tic;
  [wp, stats] = phip(TF, P62.Problem.A, ones(N, p+1), sqrt(eps), 0);
  cpu1 = toc;
  cpuphip(i) = cpu1;
end  
avecpuphip = mean(cpuphip);
fprintf('%5.2f cpu time\n',  avecpuphip);
err = norm((wexact-wp)./wexact)

% PHIPM
fprintf('------------------------------\n');
fprintf('*** PHIPM ***\n');
for i = 1:NumExp
  tic;
  [wp, stats] = phipm(TF, P62.Problem.A, ones(N, p+1), sqrt(eps), 0);
  cpu1 = toc;
  cpuphipm(i) = cpu1;
end  
fprintf('%5.2f speedup over phip\n',  avecpuphip/mean(cpuphipm));
err = norm((wexact-wp)./wexact)

fprintf('------------------------------\n');

% ===========================================================
% 6.3 - A Hermitian example
fprintf('*** 6.3 - A Hermitian example HB/bcspwr10 ***\n');
P63 = load('bcspwr10');
T0 = 0;
TF = 2;
N = 5300;
b = [ones(1,p+1); zeros(N-2, p+1); ones(1,p+1)];
b = ones(N,p+1);

% PHIPM
fprintf('------------------------------\n');
fprintf('*** Exact ***\n');
K = spdiags(ones(p, 1), 1, p, p);
A = [-P63.Problem.A, ones(N, p); zeros(p, N), K];
ep = zeros(p,1); ep(end) = 1;
bb = [ones(N, 1); ep];
[wexact2, stats] = phipm(TF, -P63.Problem.A, b, 1e-10);
[wexact,tvals,mv] = expmv_tspan(A, bb, T0, TF, 100, 'double');
wexact = wexact(1:N,end);
norm(wexact2-wexact)

% PHIP
fprintf('------------------------------\n');
fprintf('*** PHIP ***\n');
for i = 1:NumExp
  tic;
  [wp, stats] = phip(TF, -P63.Problem.A, b, 1e-5);
  cpu1 = toc;
  cpuphip(i) = cpu1;
end    
avecpuphip = mean(cpuphip);
fprintf('%5.2f cpu time\n', avecpuphip);
err = norm((wexact-wp)./wexact)

% PHIPM
fprintf('------------------------------\n');
fprintf('*** PHIPM ***\n');
for i = 1:NumExp
  tic;
  [wp, stats] = phipm(TF, -P63.Problem.A, b, 1e-5);
  cpu1 = toc;
  cpuphipm(i) = cpu1;
end      
fprintf('%5.2f sppedup over phip\n', avecpuphip/mean(cpuphipm));
err = norm((wexact-wp)./wexact)

fprintf('------------------------------\n');

% ===========================================================
% 6.4 - A forward example
fprintf('*** 6.4 - A forward example HB/gr_30_30 ***\n');
P64 = load('gr_30_30');
T0 = 0;
TF = 2;
N = 900;

% PHIPM
fprintf('------------------------------\n');
fprintf('*** Exact ***\n');
K = spdiags(ones(p, 1), 1, p, p);
A = [P64.Problem.A, ones(N, p); zeros(p, N), K];
ep = zeros(p, 1); ep(end) = 1;
b = [ones(N, 1); ep];
[wexact2, stats] = phipm(TF, P64.Problem.A, ones(N, p+1), 1e-14, 1);
[wexact,tvals,mv] = expmv_tspan(A, b, T0, TF, 100, 'double');
wexact = wexact(1:N,end);
norm(wexact2-wexact)

% PHIP
fprintf('------------------------------\n');
fprintf('*** PHIP ***\n');
for i = 1:NumExp
  tic;
  [wp, stats] = phip(TF, P64.Problem.A, ones(N, p+1), sqrt(eps), true);
  cpu1 = toc;
  cpuphip(i) = cpu1;
end
avecpuphip = mean(cpuphip)
fprintf('%5.2f cputime\n', avecpuphip);
err = norm((wexact-wp)./wexact)

% PHIPM
fprintf('------------------------------\n');
fprintf('*** PHIPM ***\n');
for i = 1:NumExp
  tic;
  [wp, stats] = phipm(TF, P64.Problem.A, ones(N, p+1), sqrt(eps), true);
  cpu1 = toc;
  cpuphipm(i) = cpu1;
end  
fprintf('%5.2f speedup over expv\n', avecpuphip/mean(cpuphipm));
err = norm((wexact-wp)./wexact)

fprintf('------------------------------\n');

% ===========================================================
% 6.7 - More terms example
fprintf('*** 6.9 - Helmoltz equation on unit square helm2d03 ***\n');
P67 = load('helm2d03');
T0 = 0;
TF = 2;
N = 392257;

% PHIPM
fprintf('------------------------------\n');
fprintf('*** Exact ***\n');
K = spdiags(ones(p, 1), 1, p, p);
A = [P67.Problem.A, ones(N, p); zeros(p, N), K];
ep = zeros(p, 1); ep(end) = 1;
b = [ones(N, 1); ep];
[wexact2, stats] = phipm(TF, P67.Problem.A, ones(N, p+1), 1e-14, true);
[wexact,tvals,mv] = expmv_tspan(A, b, T0, TF, 100, 'double');
wexact = wexact(1:N,end);
norm(wexact2-wexact)

% PHIP
fprintf('------------------------------\n');
fprintf('*** PHIP ***\n');
for i = 1:NumExp
  tic;
  [wp, stats] = phip(TF, P67.Problem.A, ones(N, p+1), sqrt(eps), true);
  cpu1 = toc;
  cpuphip(i) = cpu1;
end
avecpuphip = mean(cpuphip);
fprintf('%5.2f cpu time\n', avecpuphip);
err = norm((wexact-wp)./wexact)

% PHIPM
fprintf('------------------------------\n');
fprintf('*** PHIPM ***\n');
for i = 1:NumExp
  tic;
  [wp, stats] = phipm(TF, P67.Problem.A, ones(N, p+1), sqrt(eps), true);
  cpu1 = toc;
  cpuphipm(i) = cpu1;
end
fprintf('%5.2f speedup over phip\n', avecpuphip/mean(cpuphipm));
err = norm((wexact-wp)./wexact)

fprintf('------------------------------\n');
