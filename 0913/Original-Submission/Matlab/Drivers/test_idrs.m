%
% IDRS testproblems, 
% Uses idrs.m make_system.m idrs_ex*.m (the actual tests), results.m
% mv.m m1.m and m2.m
%
%   Martin van Gijzen
%   Version August 31, 2010
%
%   This software is distributed under the
%   ACM Software Copyright and License Agreement.
%

clear all;
close all;
clc;

% Define test problem
h = 0.1;
eps = 0.02;
beta(1) = 0/sqrt(5);
beta(2) = 1/sqrt(5);
beta(3) = 2/sqrt(5);
r = 6;
% Generate system
[A, b] = make_system( eps, beta, r, h );

n = length(b);
choice = 0;

while choice ~= 12 
   choice = menu('IDR(s) tests: ', 'No preconditioning', ...
                 'SSOR preconditioning with Eisenstat', ...
                 'SSOR preconditioning, matrices', ...
                 'SSOR preconditioning, functions', ...
                 'Residual smoothing', ...
                 'Different computations of parameter omega', ...
                 'Complex shadow vectors', ...
		 'High tolerance', ...
                 'Residual replacement', ...
		 'Restarting',...
                 'Equivalence of Bi-CGSTAB and IDR(1)', ...
		 'STOP' );
   if ( choice == 12 )
      return;
   end

   % Initialise figure:
   scrsz = get(0,'ScreenSize');
   figure('Position',[scrsz(1) + scrsz(3)/2 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);
   hold on;
   xlabel('Number of MATVECS')
   ylabel('|r|/|b|')
   grid on;
   title(['IDR(s) EXAMPLE ', num2str(choice)]);
   colour = 0;
   clear options replacements flag relres resvec M1 M2 K; 

   if ( choice == 1 )
      idrs_ex1;
   elseif ( choice == 2 )
      idrs_ex2;
   elseif ( choice == 3 )
      idrs_ex3;
   elseif ( choice == 4 )
      idrs_ex4;
   elseif ( choice == 5 )
      idrs_ex5;
   elseif ( choice == 6 )
      idrs_ex6;
   elseif ( choice == 7 )
      idrs_ex7;
   elseif ( choice == 8 )
      idrs_ex8;
   elseif ( choice == 9 )
      idrs_ex9;
   elseif ( choice == 10 )
      idrs_ex10;
   elseif ( choice == 11 )
      idrs_ex11;
   end
   disp('===============================================================================================');

   hold off;
end
