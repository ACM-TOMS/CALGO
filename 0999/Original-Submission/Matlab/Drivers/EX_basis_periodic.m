% EXAMPLE 2: periodic MDB-spline basis

% Construction of a B-spline multi-patch
P1 = B_patch(3, [0, 2]);
P2 = B_patch(4, [0, 3/2, 4], 2);
P3 = B_patch(5, [0, 3]);
MP = MDB_patch([P1, P2, P3]);
m = 901;

% Construction and visualization of C^2 periodic MDB-splines
rr = [2, 2];  rp = 3;
figure(1); clf;
H = MDB_extraction_periodic(MP, rr, rp);
MDB_visualization_all(MP, H, m, 'LineWidth', 2);
axis([0 9 0 1]);
title('C^2 periodic MDB-splines of degrees (3,4,5)');

% Construction and visualization of related derivatives
figure(2); clf;
xx = linspace(0, 9, m);
M = MDB_differentiation_all(MP, H, 1, xx);
plot(xx, M, 'LineWidth', 2);
axis([0 9 -1 1]);
title('Derivative of C^2 periodic MDB-splines of degrees (3,4,5)');
