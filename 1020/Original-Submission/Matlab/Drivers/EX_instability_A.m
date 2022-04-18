% EXAMPLE 2/A: Instability of ECT-spaces
% Null-spaces of different dimensions

% Visualization of Bernstein functions
p = [9, 10];
alpha0 = 1/(6*pi);
alpha1 = 1/(3*pi);
roots = complex([0, alpha0, alpha1, alpha0], [1, 0, 0, 1]);
xi = [11*pi/2, 49*pi/8];
m = 501;
xx = linspace(xi(1), xi(2), m);
for pind = 1:length(p)
   P = TB_patch_tcheb(p(pind), xi, roots);
   figure(pind); clf;
   TB_visualization_all(P, m, 'LineWidth', 2);
   title(num2str(p(pind), '$p = %d$'), 'Interpreter', 'latex');
   axis([xi, 0, 1]);
   M = TB_evaluation_all(P, xx);
   pouv = max(abs(sum(M, 1) - 1));
   minv = min(M(:));
   fprintf('p = %d, min = %.3e, pou = 1+%.3e\n', p(pind), minv, pouv);
end
