% EXAMPLE 1: MDB-spline basis

% Construction of a B-spline multi-patch
P1 = B_patch(3, [0, 2]);
P2 = B_patch(4, [0, 3/2, 4], 2);
P3 = B_patch(5, [0, 3]);
MP = MDB_patch([P1, P2, P3]);
m = 901;

% Construction and visualization of C^0 MDB-splines
rr = [0, 0];
figure(1); clf;
H0 = MDB_extraction(MP, rr);
MDB_visualization_all(MP, H0, m, 'LineWidth', 2);
axis([0 9 0 1]);
title('C^0 MDB-splines of degrees (3,4,5)');

% Construction and visualization of C^1 MDB-splines
rr = [1, 1];
figure(2); clf;
H1 = MDB_extraction(MP, rr);
MDB_visualization_all(MP, H1, m, 'LineWidth', 2);
axis([0 9 0 1]);
title('C^1 MDB-splines of degrees (3,4,5)');

% Construction and visualization of C^2 MDB-splines
rr = [2, 2];
figure(3); clf;
H2 = MDB_extraction(MP, rr);
MDB_visualization_all(MP, H2, m, 'LineWidth', 2);
axis([0 9 0 1]);
title('C^2 MDB-splines of degrees (3,4,5)');

% Construction and visualization of related derivatives
figure(4); clf;
xx = linspace(0, 9, m);
M = MDB_differentiation_all(MP, H0, 1, xx);
plot(xx, M, 'LineWidth', 2);
axis([0 9 -3 3]);
title('First derivative of C^0 MDB-splines of degrees (3,4,5)');
figure(5); clf;
M = MDB_differentiation_all(MP, H1, 1, xx);
plot(xx, M, 'LineWidth', 2);
axis([0 9 -3 3]);
title('First derivative of C^1 MDB-splines of degrees (3,4,5)');
figure(6); clf;
M = MDB_differentiation_all(MP, H2, 1, xx);
plot(xx, M, 'LineWidth', 2);
axis([0 9 -3 3]);
title('First derivative of C^2 MDB-splines of degrees (3,4,5)');
