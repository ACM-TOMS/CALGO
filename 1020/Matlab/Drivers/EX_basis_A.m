% EXAMPLE 1/A: MDTB-spline basis

% Construction of four different ECT-spaces
P1 = TB_patch_poly(3, [0, 1]);
P2 = TB_patch_gexp(4, [0, 1], 3);
P3 = TB_patch_gtrig(4, [0, 1], 3/2);
P4 = TB_patch_tcheb(6, [0, 1], complex([1, -1, 0], [0, 0, 2]));
PP = [P1, P2, P3, P4];
MP = MDTB_patch(PP);
m = 501;

% Construction and visualization of Bernstein basis
figure(1); clf; 
for i=1:4
   subplot(2, 2, i);
   TB_visualization_all(PP(i), m, 'LineWidth', 2);
   axis([0, 1, 0, 1]);
end
xh = axes();
th = title('Bernstein basis functions');
set(xh, 'Visible', 'off', 'HitTest', 'on');
set(th, 'Visible', 'on', 'Position', [0.5, 0.45]);

% Construction and visualization of MDTB-splines
figure(2); clf;
rr = [2, 3, 3];
H = MDTB_extraction(MP, rr);
MDTB_visualization_all(MP, H, m, 'LineWidth', 2);
title('MDTB-splines');
axis([0, 4, 0, 1]);

% Construction and visualization of periodic MDTB-splines
figure(3); clf;
rr = [2, 3, 3];  rper = 2;
Hper = MDTB_extraction_periodic(MP, rr, rper);
MDTB_visualization_all(MP, Hper, m, 'LineWidth', 2);
title('periodic MDTB-splines');
axis([0, 4, 0, 1]);
