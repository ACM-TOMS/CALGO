% EXPERIMENT1 - Produce results given in Table 1.
%
% DESCRIPTION:
%   Compare phipm, phip, expv and phiv on various large sparse matrices.

% ===========================================================
NumExp = 100;
% ===========================================================
% 6.2 - A nonsymmetric example
fprintf('*** 6.2 - A nonsymmetric example orani678 ***\n');
P62 = load('orani678')
T0 = 0;
TF = 10;
N = 2529;
b = ones(N, 1);

% EXACT
fprintf('------------------------------\n');
fprintf('*** Exact ***\n');
[wexact2, stats] = phipm(TF, P62.Problem.A, b, eps, false);
[wexact,tvals,mv] = expmv_tspan(P62.Problem.A, b, T0, TF, 100, 'double');
wexact = wexact(1:N,end);
norm(wexact2-wexact)

% EXPV
fprintf('------------------------------\n');
fprintf('*** EXPV ***\n');
for i = 1:NumExp
  cpu0 = cputime;
  [we, err] = expv(TF, P62.Problem.A, b, sqrt(eps));
  cpu1 = cputime;
  cpuexpv(i) = cpu1 - cpu0;
end
avecpuexpv = mean(cpuexpv);
fprintf('%5.2f cpu time\n', avecpuexpv);
err = norm((wexact-we)./wexact)

% PHIV
fprintf('------------------------------\n');
fprintf('*** PHIV ***\n');
for i = 1:NumExp
  cpu0 = cputime;
  [we, err] = phiv(TF, P62.Problem.A, zeros(N, 1), b, sqrt(eps));
  cpu1 = cputime;
  cpuphiv(i) = cpu1 - cpu0;
end
fprintf('%5.2f speedup over expv\n', avecpuexpv/mean(cpuphiv));
err = norm((wexact-we)./wexact)

% PHIP
fprintf('------------------------------\n');
fprintf('*** PHIP ***\n');
for i = 1:NumExp
  cpu0 = cputime;
  [wp, stats] = phip(TF, P62.Problem.A, b, sqrt(eps), false);
  cpu1 = cputime;
  cpuphip(i) = cpu1 - cpu0;
end  
fprintf('%5.2f speedup over expv\n',  avecpuexpv/mean(cpuphip));
err = norm((wexact-wp)./wexact)
stats

% PHIPM
fprintf('------------------------------\n');
fprintf('*** PHIPM ***\n');
for i = 1:NumExp
  cpu0 = cputime;
  [wp, stats] = phipm(TF, P62.Problem.A, b, sqrt(eps), false);
  cpu1 = cputime;
  cpuphipm(i) = cpu1 - cpu0;
end  
fprintf('%5.2f speedup over expv\n',  avecpuexpv/mean(cpuphipm));
err = norm((wexact-wp)./wexact)

fprintf('------------------------------\n');

% ===========================================================
% 6.3 - A Hermitian example
fprintf('*** 6.3 - A Hermitian example HB/bcspwr10 ***\n');
P63 = load('bcspwr10');
T0 = 0;
TF = 2;
N = 5300;
b = [1; zeros(N-2, 1); 1];

% PHIPM
fprintf('------------------------------\n');
fprintf('*** Exact ***\n');
[wexact2, stats] = phipm(TF, -P63.Problem.A, b, 1e-10, true);
[wexact,tvals,mv] = expmv_tspan(-P63.Problem.A, b, T0, TF, 100, 'double');
wexact = wexact(1:N,end);
I = find(abs(wexact) > 1e-5);
norm(wexact2-wexact)

% EXPV
fprintf('------------------------------\n');
fprintf('*** EXPV ***\n');
for i = 1:NumExp
  cpu0 = cputime;
  [we, err] = expv(TF, -P63.Problem.A, b, 1e-5);
  cpu1 = cputime;
  cpuexpv(i) = cpu1 - cpu0;
end  
avecpuexpv = mean(cpuexpv);
fprintf('%5.2f cpu time\n', avecpuexpv);
err = norm((wexact(I)-we(I))./wexact(I))

% PHIV
fprintf('------------------------------\n');
fprintf('*** PHIV ***\n');
for i = 1:NumExp
  cpu0 = cputime;
  [we, err] = phiv(TF, -P63.Problem.A, zeros(N,1), b, 1e-5);
  cpu1 = cputime;
  cpuphiv(i) = cpu1 - cpu0;
end  
fprintf('%5.2f speedup over expv\n', avecpuexpv/mean(cpuphiv));
err = norm((wexact(I)-we(I))./wexact(I))

% PHIP
fprintf('------------------------------\n');
fprintf('*** PHIP ***\n');
for i = 1:NumExp
  cpu0 = cputime;
  [wp0, stats] = phip(TF, -P63.Problem.A, b, 1e-5, true);
  cpu1 = cputime;
  cpuphip(i) = cpu1 - cpu0;
end  
fprintf('%5.2f speedup over expv\n', avecpuexpv/mean(cpuphip));
err = norm((wexact(I)-wp0(I))./wexact(I))

% PHIPM
fprintf('------------------------------\n');
fprintf('*** PHIPM ***\n');
for i = 1:NumExp
  cpu0 = cputime;
  [wp, stats] = phipm(TF, -P63.Problem.A, b, 1e-5, true);
  cpu1 = cputime;
  cpuphipm(i) = cpu1 - cpu0;
end      
fprintf('%5.2f speedup over expv\n', avecpuexpv/mean(cpuphipm));
err = norm((wexact(I)-wp(I))./wexact(I))

fprintf('------------------------------\n');

% ===========================================================
% 6.4 - A forward-backward example
fprintf('*** 6.4 - A forward-backward example HB/gr_30_30 ***\n');
P64 = load('gr_30_30');
TF = 2;
N = 900;
wexact = ones(N, 1);

% EXPV
fprintf('------------------------------\n');
fprintf('*** EXPV ***\n');
for i = 1:NumExp
  cpu0 = cputime;
  [we, err] = expv(TF, P64.Problem.A, wexact, 1e-14);
  [we, err] = expv(-TF, P64.Problem.A, we, 1e-14);
  cpu1 = cputime;
  cpuexpv(i) = cpu1 - cpu0;
end
avecpuexpv = mean(cpuexpv);
fprintf('%5.2f cpu time\n', avecpuexpv);
err = norm((wexact-we)./wexact)

% PHIV
fprintf('------------------------------\n');
fprintf('*** PHIV ***\n');
for i = 1:NumExp
  cpu0 = cputime;
  [we, err] = phiv(TF, P64.Problem.A, zeros(N, 1), wexact, 1e-14);
  [we, err] = phiv(-TF, P64.Problem.A, zeros(N, 1), we, 1e-14);
  cpu1 = cputime;
  cpuphiv(i) = cpu1 - cpu0;
end
fprintf('%5.2f speedup over expv\n', avecpuexpv/mean(cpuphiv));
err = norm((wexact-we)./wexact)

% PHIP
fprintf('------------------------------\n');
fprintf('*** PHIP ***\n');
for i = 1:NumExp
  cpu0 = cputime;
  [wp, stats] = phip(TF, P64.Problem.A, wexact, 1e-14, true);
  [wp, stats] = phip(-TF, P64.Problem.A, wp, 1e-14, true);
  cpu1 = cputime;
  cpuphip(i) = cpu1 - cpu0;
end
fprintf('%5.2f speedup over expv\n', avecpuexpv/mean(cpuphip));
err = norm((wexact-wp)./wexact)
stats

% PHIPM
fprintf('------------------------------\n');
fprintf('*** PHIPM ***\n');
for i = 1:NumExp
  cpu0 = cputime;
  [wp, stats] = phipm(TF, P64.Problem.A, wexact, 1e-14, true);
  [wp, stats] = phipm(-TF, P64.Problem.A, wp, 1e-14, true);
  cpu1 = cputime;
  cpuphipm(i) = cpu1 - cpu0;
end  
fprintf('%5.2f speedup over expv\n', avecpuexpv/mean(cpuphipm));
err = norm((wexact-wp)./wexact)

fprintf('------------------------------\n');

% ===========================================================
% 6.7 - More terms example
fprintf('*** 6.7 - A non-homogenous example helm2d03 ***\n');
P67 = load('helm2d03');
T0 = 0;
TF = 2;
p = 1;
N = 392257;

% PHIPM
fprintf('------------------------------\n');
fprintf('*** Exact ***\n');
K = spdiags(ones(p, 1), 1, p, p);
A = [P67.Problem.A, ones(N, p); zeros(p, N), K];
ep = zeros(p, 1); ep(end) = 1;
b = [ones(N, 1); ep];
[wexact2, stats] = phipm(TF, P67.Problem.A, ones(N, 2), 1e-14, true);
[wexact,tvals,mv] = expmv_tspan(A, b, T0, TF, 100, 'double');
wexact = wexact(1:N,end);
norm(wexact2-wexact)

% PHIV
fprintf('------------------------------\n');
fprintf('*** PHIV ***\n');
for i = 1:NumExp
  cpu0 = cputime;
  [we, err] = phiv(TF, P67.Problem.A, ones(N, 1), ones(N, 1), sqrt(eps));
  cpu1 = cputime;
  cpuphiv(i) = cpu1 - cpu0;
end
avecpuphiv = mean(cpuphiv);
fprintf('%5.2f cpu time\n', avecpuphiv);
err = norm((wexact-we)./wexact)

% PHIP
fprintf('------------------------------\n');
fprintf('*** PHIP ***\n');
for i = 1:NumExp
  cpu0 = cputime;
  [wp, stats] = phip(TF, P67.Problem.A, ones(N, 2), sqrt(eps), true);
  cpu1 = cputime;
  cpuphip(i) = cpu1 - cpu0;
end
fprintf('%5.2f speedup over phiv\n', avecpuphiv/mean(cpuphip));
err = norm((wexact-wp)./wexact)

% PHIPM
fprintf('------------------------------\n');
fprintf('*** PHIPM ***\n');
for i = 1:NumExp
  cpu0 = cputime;
  [wp, stats] = phipm(TF, P67.Problem.A, ones(N, 2), sqrt(eps), true);
  cpu1 = cputime;
  cpuphipm(i) = cpu1 - cpu0;
end
fprintf('%5.2f speedup over phiv\n', avecpuphiv/mean(cpuphipm));
err = norm((wexact-wp)./wexact)

fprintf('------------------------------\n');
