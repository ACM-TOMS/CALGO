% EXAMPLE 2/B: Instability of ECT-spaces
% Generalized polynomial spaces with trigonometric functions

% Visualization of Bernstein functions
w = 1/3;
p = 10;
m = 501;
xx = linspace(0, 1, m);

figure(1); clf;
P1 = TB_patch_gtrig(p, [0, 1], w);
TB_visualization_all(P1, m, 'LineWidth', 2);
title('TB\_patch\_gtrig');
axis([0, 1, 0, 1]);
M1 = TB_evaluation_all(P1, xx);
minv1 = min(M1(:));
pouv1 = max(abs(sum(M1, 1) - 1));
fprintf('gtrig: min = %.3e, pou = 1+%.3e\n', minv1, pouv1);

figure(2); clf;
P2 = TB_patch_tcheb(p, [0, 1], complex(0, w));
TB_visualization_all(P2, m, 'LineWidth', 2);
title('TB\_patch\_tcheb');
axis([0, 1, 0, 1]);
M2 = TB_evaluation_all(P2, xx);
minv2 = min(M2(:));
pouv2 = max(abs(sum(M2, 1) - 1));
fprintf('tcheb: min = %.3e, pou = 1+%.3e\n', minv2, pouv2);
