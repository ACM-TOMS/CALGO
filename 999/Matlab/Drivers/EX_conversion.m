% EXAMPLE 3: MDB-spline conversion

% Construction of a B-spline multi-patch (2x)
P1 = B_patch(7, [0, 1]);
P2 = B_patch(2, [0, 1]);
P3 = B_patch(3, [0, 1]);
MP = MDB_patch([P1, P2, P3]);
MPfull = MDB_patch([P1, P1, P1]);
m = 301;

% Construction and visualization of MDB-splines
rr = [2, 1];
figure(1); clf;
H = MDB_extraction(MP, rr);
MDB_visualization_all(MP, H, m, 'LineWidth', 2);
title('MDB-splines of degrees (7,2,3)');

% Construction and visualization of B-splines (2x)
figure(2); clf;
Hfull = MDB_extraction(MPfull, rr);
MDB_visualization_all(MPfull, Hfull, m, 'LineWidth', 2);
title('(MD)B-splines of degree 7');
figure(3); clf;
Pfull = B_patch(7, [0, 1, 2, 3], rr);
B_visualization_all(Pfull, m, 'LineWidth', 2);
title('B-splines of degree 7');

% Conversion and visualization of a C^1 spline function (3x)
cc = [7, 4, 10, 1, 4, 5/2, 2, 3/2, 2, 3];
ccfull = MDB_conversion(MPfull, Hfull, MP, H, cc);
figure(4); clf;
MDB_visualization_spline(MP, H, cc, m, 'LineWidth', 2);
title('C^1 spline function in MDB-spline representation');
figure(5); clf;
MDB_visualization_spline(MPfull, Hfull, ccfull, m, 'LineWidth', 2);
title('C^1 spline function in (MD)B-spline representation');
figure(6); clf;
B_visualization_spline(Pfull, ccfull, m, 'LineWidth', 2);
title('C^1 spline function in B-spline representation');
