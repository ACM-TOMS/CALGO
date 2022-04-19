% EXAMPLE 3/A: Critical length of ECT-spaces
% Generalized polynomial spaces with trigonometric functions

% Computation of guesses of critical lengths
w = 1;
p = 2:10;
l = 2:0.001:20;
xx = 0:0.001:l(end);
mp = zeros(length(p), length(l));
lp = zeros(length(p), 1);
tol = 1e-6;
for pind = 1:length(p)
   for lind = 1:length(l)
      P = TB_patch_gtrig(p(pind), [0, l(lind)], w);
      M = TB_evaluation_all(P, xx);
      minv = min(M(:));
      if minv < -tol
         lm = lind - 1;
         lp(pind) = l(lm);
         fprintf('l''_%d = %.3f\n', p(pind),  lp(pind));
         break;
      else
         mp(pind, lind) = minv;
      end
   end
end

% Visualization of parameterizations of critical lengths
figure(1); clf;
m = 0:0.01:10;
plot(m, lp * (1./m), 'LineWidth', 2);
title('guess of critical length');
xlabel('$\beta$', 'Interpreter', 'latex');
ylabel('$\ell''_p$', 'Interpreter', 'latex');
legend(cellstr(num2str(p', '$p = %d$')), 'Interpreter', 'latex');
axis([0, 10, 0, 10]);
