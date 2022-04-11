% EXAMPLE 1/B: MDTB-spline basis
% Representation of a square with circular corners

% Construction of two different ECT-spaces
s = 0;
l = 4^s;
P1 = TB_patch_ptrig(2, [0, pi/2], 1);
P2 = TB_patch_poly(1, [0, l]);
PP = [P1, P2, P1, P2, P1, P2, P1, P2];
MP = MDTB_patch(PP);
m = 501;

% Construction and visualization of C^1 MDTB-splines
figure(1); clf;
Hper = MDTB_extraction_periodic(MP, 1, 1);
MDTB_visualization_all(MP, Hper, m, 'LineWidth', 2);
title('periodic MDTB-splines');
axis([0, 4*l+2*pi, 0, 1]);

% Construction and visualization of C^1 MDTB-spline curve
figure(2); clf;
cc = [1, -1, -1, 1; 1, 1, -1, -1];
MDTB_visualization_curve(MP, Hper, cc, m, 'LineWidth', 2);
patch(cc(1,:), cc(2,:), 0, 'FaceColor', 'none', 'Marker', 'p');
title('square with circular corners');
axis([-1.3, 1.3, -1.3, 1.3]);
axis equal;
