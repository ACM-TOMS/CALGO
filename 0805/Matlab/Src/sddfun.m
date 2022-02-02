function sddfun(A)

% This is a test to compare the SDD with the SVD and to compare different
% SDD starting strategies.  It computes the rank of the matrix A and its
% SVD. It then computes a k-term SDD with each of four different starting
% strategies where k is the rank of the matrix. 
%
% A table and four plots provide comparisons. The table compares the
% different starting strategies with respect to final relative residual,
% average number of inner iterations, and sparsity of factors. The plots
% compare residual norm vs. iteration and  residual norm vs. storage for the
% SVD and the 4 SDDs. We also compare the singular values with the (scaled)
% diagonal value from the SDDs. Lastly, we compare work vs. the residual
% norm for each SDD method.
%
%
%SDDPACK: Software for the Semidiscrete Decomposition.
%Copyright (c) 1999 Tamara G. Kolda and Dianne P. O'Leary. 

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the Free 
% Software Foundation; either version 2 of the License, or (at your option)
% any later version.  
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
% or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
% for more details.  
%
% You should have received a copy of the GNU General Public License along
% with this program; if not, write to the Free Software Foundation, Inc., 59
% Temple Place - Suite 330, Boston, MA 02111-1307, USA.  


%----------------------------------------------------------------------
% Initialization
%----------------------------------------------------------------------
[m, n] = size(A);			% size of A
rho0 = sum(sum(A .* A));		% 2-norm of A

%----------------------------------------------------------------------
% SVD Stuff
%----------------------------------------------------------------------
fprintf(1, '\n');
disp('Computing SVD');
[U, S, V] = svd(full(A));		% svd of A
S = diag(S);				% diagonalize S
k = size(S,1);				% rank of A

fprintf(1, '\n');
fprintf(1, 'Rank of matrix: %d\n', k);

clear svdrho				% residual or SVD
svdrho(k) = 0;				% final residual
for i = k-1 : -1 : 1			% working backwards...
  svdrho(i) = svdrho(i+1) + S(i+1)^2;
end

svdstor = 4 * (m + n + 1) * [1 : k];	% storage for SVD

%----------------------------------------------------------------------
% SDD Stuff
%----------------------------------------------------------------------

fprintf(1, '\n');
disp('Computing SDD with option 1');
[D1, X1, Y1, its1, rho1, iits1] = sdd(A, k, 0.01, 100, 0.0, 1);
disp('Computing SDD with option 2');
[D2, X2, Y2, its2, rho2] = sdd(A, k, 0.01, 100, 0.0, 2);
disp('Computing SDD with option 3');
[D3, X3, Y3, its3, rho3] = sdd(A, k, 0.01, 100, 0.0, 3);
disp('Computing SDD with option 4');
[D4, X4, Y4, its4, rho4] = sdd(A, k, 0.01, 100, 0.0, 4);

its1 = its1 + (iits1 - 1) / 2.0;

clear Dplus1
for i = 1 : k
  Dplus1(i) = D1(i) * norm(X1(:,i)) * norm(Y1(:,i));
end
clear Dplus2
for i = 1 : k
  Dplus2(i) = D2(i) * norm(X2(:,i)) * norm(Y2(:,i));
end
clear Dplus3
for i = 1 : k
  Dplus3(i) = D3(i) * norm(X3(:,i)) * norm(Y3(:,i));
end
clear Dplus4
for i = 1 : k
  Dplus4(i) = D4(i) * norm(X4(:,i)) * norm(Y4(:,i));
end

sddstor = (4 + .25 * (m + n)) * [1 : k]; % storage for SDD

density1 = 100 * (sum(sum(abs(X1))) + sum(sum(abs(Y1)))) / (2 * m * n);
density2 = 100 * (sum(sum(abs(X2))) + sum(sum(abs(Y2)))) / (2 * m * n);
density3 = 100 * (sum(sum(abs(X3))) + sum(sum(abs(Y3)))) / (2 * m * n);
density4 = 100 * (sum(sum(abs(X4))) + sum(sum(abs(Y4)))) / (2 * m * n);

clear sits1 sits2 sits3 sits4
sits1(1) = its1(1);
sits2(1) = its2(1);
sits3(1) = its3(1);
sits4(1) = its4(1);

for i = 2 : k
  sits1(i) = sits1(i-1) + its1(i);
  sits2(i) = sits2(i-1) + its2(i);
  sits3(i) = sits3(i-1) + its3(i);
  sits4(i) = sits4(i-1) + its4(i);
end

%----------------------------------------------------------------------
% Plotting
%----------------------------------------------------------------------

c = ['b' 'm' 'r' 'c' 'g'];		% colors
l = ['- '; ': '; '-.'; '--'; '- '];	% line pattern
d = ['o'; '^'; '+'; '*'; 'x'];		% dot pattern

fprintf(1,'\n');
fprintf(1, 'Method Color Line Dot\n');
fprintf(1, '------ ----- ---- ---\n');
fprintf(1, 'SVD      %1s    %2s   %1s\n', c(5), l(5,:), d(5));
for i = 1:4
  fprintf(1, 'SDD-%d    %1s    %2s   %1s\n', i, c(i), l(i,:), d(i));
end
fprintf(1,'\n');

%%% Residual Norms vs. No. of Terms
subplot(2,2,1);
plot([0:k], sqrt([rho0 svdrho]/rho0), [c(5) l(5,:)]);
hold on;
plot([0:k], sqrt([rho0 rho1]/rho0), [c(1) l(1,:)]);
plot([0:k], sqrt([rho0 rho2]/rho0), [c(2) l(2,:)]);
plot([0:k], sqrt([rho0 rho3]/rho0), [c(3) l(3,:)]);
plot([0:k], sqrt([rho0 rho4]/rho0), [c(4) l(4,:)]);
hold off;
title('Residual Norms vs. No. of Terms');

%%% Singular vs. SDD values
subplot(2,2,2);
plot([1:k], sqrt(S), [c(5) d(5)]);
hold on;
plot([1:k], sqrt(Dplus1), [c(1) d(1)]);
plot([1:k], sqrt(Dplus2), [c(2) d(2)]);
plot([1:k], sqrt(Dplus3), [c(3) d(3)]);
plot([1:k], sqrt(Dplus4), [c(4) d(4)]);
hold off;
title('Singular vs. SDD values');

%%% Residual Norm vs. Storage
subplot(2, 2, 3);
semilogx(svdstor, sqrt(svdrho/rho0), [c(5) l(5,:)]);
hold on;
semilogx(sddstor, sqrt(rho1/rho0), [c(1) l(1,:)]);
semilogx(sddstor, sqrt(rho2/rho0), [c(2) l(2,:)]);
semilogx(sddstor, sqrt(rho3/rho0), [c(3) l(3,:)]);
semilogx(sddstor, sqrt(rho4/rho0), [c(4) l(4,:)]);
hold off;
title('Residual Norms vs. Storage');

%%% SDD Inner Iterations
subplot(2, 2, 4);
plot(sits1, sqrt([rho1]/rho0), [c(1) l(1,:)]);
hold on;
plot(sits2, sqrt([rho2]/rho0), [c(2) l(2,:)]);
plot(sits3, sqrt([rho3]/rho0), [c(3) l(3,:)]);
plot(sits4, sqrt([rho4]/rho0), [c(4) l(4,:)]);
hold off;
title('Residual Norms vs. No. Inner Iterations');

fprintf(1, 'Method  Rel Resid  Inn Its  Density\n');
fprintf(1, '------  ---------  -------  -------\n');
fprintf(1, 'SDD-1     %5.2f     %5.2f    %5.2f\n', 100 * sqrt(rho1(k)/rho0), mean(its1), density1);
fprintf(1, 'SDD-2     %5.2f     %5.2f    %5.2f\n', 100 * sqrt(rho2(k)/rho0), mean(its2), density2);
fprintf(1, 'SDD-3     %5.2f     %5.2f    %5.2f\n', 100 * sqrt(rho3(k)/rho0), mean(its3), density3);
fprintf(1, 'SDD-4     %5.2f     %5.2f    %5.2f\n', 100 * sqrt(rho4(k)/rho0), mean(its4), density4);
