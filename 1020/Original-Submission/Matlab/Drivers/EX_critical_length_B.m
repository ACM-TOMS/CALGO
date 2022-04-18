% EXAMPLE 3/B: Critical length of ECT-spaces
% Null-spaces with only trigonometric functions

% Computation of guesses of critical lengths
p = 6;
w = 3:0.1:10;
l = 0.1:0.001:5;
xx = 0:0.001:l(end);
mp = zeros(length(w), length(l));
lp = zeros(length(w), 1);
skip = 10;
tol = 1e-6;
for wind = 1:length(w)
   for lind = 1:length(l)
      P = TB_patch_tcheb(p, [0, l(lind)], complex(0, [1, 2, w(wind)]));
      M = TB_evaluation_all(P, xx);
      minv = min(M(:));
      if minv < -tol
         lm = lind - 1;
         lp(wind) = l(lm);
         if mod(wind, skip) == 1
            fprintf('w = %d, l''_%d = %.3f\n', w(wind), p,  lp(wind));
         end
         break;
      else
         mp(wind, lind) = minv;
      end
   end
end

% Visualization of guesses of critical lengths
figure(1); clf;
plot(w, lp, 'LineWidth', 2);
title('guess of critical length');
xlabel('$\beta$', 'Interpreter', 'latex');
ylabel(num2str(p, '$\\ell''_%d$'), 'Interpreter', 'latex');
