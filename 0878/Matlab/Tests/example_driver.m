% EXAMPLE_DRIVER
%
%   EXAMPLE_DRIVER demonstrates simple calls to the main functions of the
%   package, var_ll, varma_llc, varma_llm, and varma_sim. More thorough (and 
%   perhaps realistic) demonstration is provided for by the program demorun,
%   where actual parameter estimation is carried out.
%
%   The following output is obtained:
%     f1 = -32.76150
%     f2 = -29.98882
%     f3 = -33.25909
%     f4 = -30.62788
% 
%     d1 =  -0.71  0.33 -0.22 -0.25 -2.33  2.29 -3.25
%     d2 =  -1.00  0.26 -0.10 -0.56 -1.93  1.87 -3.18 -0.05  0.76
%     d3 =  -1.04 -0.14 -0.44 -0.99 -0.88 -0.42 -0.51 -1.23 -2.32  2.32 -3.23
%     d4 =  -1.31 -0.38 -0.16 -1.24 -1.23 -0.59 -0.33 -1.32 -1.94  1.91 -3.13 -0.06  0.46
% 
%     X1 = -0.959 -0.408  0.723  0.385  1.453
%          -0.261  1.087  1.354  0.610 -0.435
%     X2 = -3.976 -3.077  2.141  1.553  1.076
%          -1.605 -0.741 -0.977  1.022  2.929

function example_driver

  % DEFINE DATA MATRICES
  A = [ 0.4  0.2 ...
    ;   0.0  0.5 ];
  
  B = [ 0.3  0.1 ...
    ;   0.2  0.3 ];
  
  Sig = [ 3  1 ...
    ;     1  2 ];
  
  mu = [0.1  0.1];
  
  X = [ .5  .3  .2  .6  .1  .4  .3  .4  .7  .0  .0  .6 ...
    ;   .8  .5  .1  .8  .2  .8  .9  .2  .0  .5  .2  .7 ];
  
  miss = false(size(X));  
  miss(1,1) = true;
  miss(1,3) = true;
  
  % LIKELIHOOD CALCULATION (WITH DERIVATIVES)
  [f1, ok, d1] = var_ll(X, A, Sig);
  [f2, ok, d2] = var_ll(X, A, Sig, mu, miss);
  [f3, ok, d3] = varma_llc(X, A, B, Sig);
  [f4, ok, d4] = varma_llm(X, A, B, Sig, mu, miss);
  %
  fprintf('f1 = %.5f\nf2 = %.5f\nf3 = %.5f\nf4 = %.5f\n', f1, f2, f3, f4)
  fprintf('\nd1 = '), fprintf(' %5.2f', d1),
  fprintf('\nd2 = '), fprintf(' %5.2f', d2),
  fprintf('\nd3 = '), fprintf(' %5.2f', d3),
  fprintf('\nd4 = '), fprintf(' %5.2f', d4), fprintf('\n\n')
  
  % SIMULATION
  randn('state', 0)
  X = varma_sim(A, B, Sig, 5, mu, 2);
  disp([['X1 = ';'     '] num2str(X(:,:,1),' %6.3f')])
  disp([['X2 = ';'     '] num2str(X(:,:,2),' %6.3f')])

end