% EXAMPLE 2/C: Instability of MDTB-spline spaces
% MDTB-spline spaces with generalized polynomial segments

% Visualization of MDTB-spline functions
xi = [-100, -99, 0, 99, 100];
tp = [2, 1, 1, 2];
w = [1/3 1/33 1/33 1/3];
m = 1001;
xx = linspace(-100, 100, m);

figure(1); clf;
pp1 = [10, 9, 9, 10];
rr1 = [8, 8, 8];
MP1 = MDTB_patch_gpoly(pp1, xi, tp, w);
H1 = MDTB_extraction(MP1, rr1);
MDTB_visualization_all(MP1, H1, m, 'LineWidth', 2);
title(['MDTB-splines of high degrees [' num2str(pp1) ']']);
axis([-100, 100, 0, 1]);
M1 = MDTB_evaluation_all(MP1, H1, xx);
M1rev = rot90(M1, 2);
minv1 = min(M1(:));
pouv1 = max(abs(sum(M1, 1) - 1));
drev1 = max(abs(M1(:) - M1rev(:)));
fprintf('mdtb high:  min = %.3e, pou = 1+%.3e\n', minv1, pouv1);
fprintf('mdtb high:  diff rev = %.3e\n', drev1);

figure(2); clf;
pp2 = [6, 5, 5, 6];
rr2 = [4, 4, 4];
MP2 = MDTB_patch_gpoly(pp2, xi, tp, w);
H2 = MDTB_extraction(MP2, rr2);
MDTB_visualization_all(MP2, H2, m, 'LineWidth', 2);
title(['MDTB-splines of lower degrees [' num2str(pp2) ']']);
axis([-100, 100, 0, 1]);
M2 = MDTB_evaluation_all(MP2, H2, xx);
M2rev = rot90(M2, 2);
minv2 = min(M2(:));
pouv2 = max(abs(sum(M2, 1) - 1));
drev2 = max(abs(M2(:) - M2rev(:)));
fprintf('mdtb lower: min = %.3e, pou = 1+%.3e\n', minv2, pouv2);
fprintf('mdtb lower: diff rev = %.3e\n', drev2);
