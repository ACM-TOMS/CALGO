% EXAMPLE 3/C: Critical length of ECT-spaces
% MDTB-spline spaces with only trigonometric functions

% Computation of guesses of critical lengths
p = 6;
m = 3;
r = 0:p-1;
w = 3:0.1:10;
l = 0.1:0.001:5;
xx = 0:0.001:m*l(end);
mp = zeros(length(r), length(w), length(l));
lp = zeros(length(r), length(w));
skip = 10;
tol = 1e-6;
for rind = 1:length(r)
   rr = r(rind) * ones(1, m-1);
   for wind = 1:length(w)
      for lind = 1:length(l)
         xi = 0:l(lind):m*l(lind);
         MP = MDTB_patch_tcheb(p, xi, {complex(0, [1, 2, w(wind)])});
         H = MDTB_extraction(MP, rr);
         M = MDTB_evaluation_all(MP, H, xx);
         minv = min(M(:));
         if minv < -tol
            lm = lind - 1;
            lp(rind, wind) = l(lm);
            if mod(wind, skip) == 1
               fprintf('r = %d, w = %d, l''_%d = %.3f\n', ...
                  r(rind), w(wind), p,  lp(rind, wind));
            end
            break;
         else
            mp(rind, wind, lind) = minv;
         end
      end
   end
end

% Visualization of guesses of critical lengths
figure(1); clf;
plot(w, lp, 'LineWidth', 2);
title('guess of critical length');
xlabel('$\beta$', 'Interpreter', 'latex');
ylabel(num2str(p, '$\\ell''_%d$'), 'Interpreter', 'latex');
legend(cellstr(num2str(r', '$r = %d$')), 'Interpreter', 'latex');
